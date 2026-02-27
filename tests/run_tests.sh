#!/usr/bin/env bash
# =============================================================================
# run_tests.sh
# ------------
# Run all tests for the augustus_bankit_toolkit.
#
# Usage:
#   bash tests/run_tests.sh          # Run all tests
#   bash tests/run_tests.sh -v       # Verbose output
#   bash tests/run_tests.sh --check  # Also run end-to-end check
#
# Requirements:
#   Python 3.6+ in PATH
#   Training data in expected location (for integration tests)
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "$SCRIPT_DIR")"
TRAINING_DIR="/prj/pflaphy-pacbio/Claud_AI/HMA4/Augustus_ to Feature_file_to_genbank_file"

# ── Colours ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'  # No colour

info()    { echo -e "${BLUE}[INFO]${NC}  $*"; }
success() { echo -e "${GREEN}[PASS]${NC}  $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
fail()    { echo -e "${RED}[FAIL]${NC}  $*"; }

VERBOSE=''
RUN_E2E=false
for arg in "$@"; do
    case "$arg" in
        -v|--verbose) VERBOSE='-v' ;;
        --check)      RUN_E2E=true ;;
    esac
done

echo "========================================================================"
echo " AUGUSTUS BankIt Toolkit – Test Suite"
echo "========================================================================"
echo " Toolkit:  $TOOLKIT_DIR"
echo " Python:   $(python3 --version 2>&1)"
echo " Date:     $(date)"
echo "========================================================================"
echo ""

PASS_COUNT=0
FAIL_COUNT=0

run_test() {
    local name="$1"
    local cmd="$2"
    echo -n "  Testing: $name ... "
    if eval "$cmd" > /tmp/_test_out.txt 2>&1; then
        success "OK"
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        fail "FAILED"
        cat /tmp/_test_out.txt
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
}

# ── Unit tests ────────────────────────────────────────────────────────────────
info "Running unit tests..."
echo ""

cd "$TOOLKIT_DIR"

run_test "AUGUSTUS parser" \
    "python3 -m pytest $VERBOSE tests/test_augustus_parser.py 2>&1"

run_test "GenBank parser" \
    "python3 -m pytest $VERBOSE tests/test_genbank_parser.py 2>&1"

run_test "GFF/GFF3/GTF parser" \
    "python3 -m pytest $VERBOSE tests/test_gff_parser.py 2>&1"

run_test "Feature table writer" \
    "python3 -m pytest $VERBOSE tests/test_feature_writer.py 2>&1"

echo ""

# ── Import check ──────────────────────────────────────────────────────────────
info "Checking module imports..."
run_test "lib imports" \
    "python3 -c 'from lib.models import GeneModel; from lib.fasta_utils import translate_cds; from lib.augustus_parser import parse_augustus_output; from lib.genbank_parser import parse_genbank; from lib.gff_parser import parse_gff; from lib.feature_writer import write_feature_table; from lib.ncbi_validator import run_all_validations; print(\"All imports OK\")'"

echo ""

# ── End-to-end tests with real data ───────────────────────────────────────────
if [ "$RUN_E2E" = true ]; then
    info "Running end-to-end tests with training data..."
    echo ""

    WALL10_HTML="$TRAINING_DIR/Augustus_ Wall10123.html"
    WALL10_FASTA="$TRAINING_DIR/Wall10.HMA4-123.fasta"
    BAC_GB="$TRAINING_DIR/Training_data/BAC_7C17_training_data_from_NCBI.gb"
    BAC_GFF3="$TRAINING_DIR/Training_data/BAC_7C_17.gff3"
    BAC_FASTA="$TRAINING_DIR/Training_data/BAC_7C17_complete_fasta.fasta"
    BAC_AUG_HTML="$TRAINING_DIR/Training_data/Augustus_ result display.html"
    OUT_DIR="/tmp/bankit_test_output"

    mkdir -p "$OUT_DIR"

    # Test 1: AUGUSTUS HTML → BankIt (Wall10)
    if [ -f "$WALL10_HTML" ] && [ -f "$WALL10_FASTA" ]; then
        run_test "augustus2bankit (Wall10 HMA4)" \
            "python3 '$TOOLKIT_DIR/augustus2bankit.py' \
                --augustus '$WALL10_HTML' \
                --fasta    '$WALL10_FASTA' \
                --organism 'Arabidopsis halleri subsp. halleri' \
                --isolate  'Wall_10_v1.0.0' \
                --product  'Zn/Cd P(IB)-type ATPase' \
                --outdir   '$OUT_DIR/wall10' \
                --prefix   wall10_hma4 \
                --primary-only 2>&1"

        if [ -f "$OUT_DIR/wall10/wall10_hma4_feature_table.tbl" ]; then
            info "Wall10 feature table:"
            head -30 "$OUT_DIR/wall10/wall10_hma4_feature_table.tbl"
        fi
    else
        warn "Wall10 test data not found – skipping Wall10 end-to-end test."
    fi

    # Test 2: GenBank → BankIt (BAC_7C17)
    if [ -f "$BAC_GB" ]; then
        run_test "genbank2bankit (BAC_7C17 GenBank)" \
            "python3 '$TOOLKIT_DIR/genbank2bankit.py' \
                --input    '$BAC_GB' \
                --format   genbank \
                --organism 'Arabidopsis halleri subsp. halleri' \
                --clone    'BAC_7C17' \
                --gene-list 'AhHMA4-1,AhHMA4-2' \
                --outdir   '$OUT_DIR/bac7c17_gb' \
                --prefix   bac7c17 2>&1"

        if [ -f "$OUT_DIR/bac7c17_gb/bac7c17_feature_table.tbl" ]; then
            info "BAC_7C17 GenBank → BankIt feature table (first 30 lines):"
            head -30 "$OUT_DIR/bac7c17_gb/bac7c17_feature_table.tbl"
        fi
    else
        warn "BAC_7C17 GenBank file not found – skipping GenBank e2e test."
    fi

    # Test 3: GFF3 → BankIt (BAC_7C17)
    if [ -f "$BAC_GFF3" ] && [ -f "$BAC_FASTA" ]; then
        run_test "genbank2bankit (BAC_7C17 GFF3)" \
            "python3 '$TOOLKIT_DIR/genbank2bankit.py' \
                --input    '$BAC_GFF3' \
                --fasta    '$BAC_FASTA' \
                --format   gff3 \
                --organism 'Arabidopsis halleri subsp. halleri' \
                --clone    'BAC_7C17' \
                --gene-list 'AhHMA4-1,AhHMA4-2' \
                --outdir   '$OUT_DIR/bac7c17_gff3' \
                --prefix   bac7c17_gff3 2>&1"
    else
        warn "BAC_7C17 GFF3/FASTA not found – skipping GFF3 e2e test."
    fi

    # Test 4: AUGUSTUS HTML (BAC_7C17) → BankIt
    if [ -f "$BAC_AUG_HTML" ] && [ -f "$BAC_FASTA" ]; then
        run_test "augustus2bankit (BAC_7C17 training HTML)" \
            "python3 '$TOOLKIT_DIR/augustus2bankit.py' \
                --augustus '$BAC_AUG_HTML' \
                --fasta    '$BAC_FASTA' \
                --organism 'Arabidopsis halleri subsp. halleri' \
                --clone    'BAC_7C17' \
                --outdir   '$OUT_DIR/bac7c17_aug' \
                --prefix   bac7c17_aug \
                --primary-only 2>&1"
    else
        warn "BAC_7C17 AUGUSTUS HTML/FASTA not found – skipping."
    fi

    # Test 5: Validate output files
    if [ -f "$OUT_DIR/bac7c17_gb/bac7c17_feature_table.tbl" ]; then
        run_test "validate_bankit (BAC_7C17 GenBank output)" \
            "python3 '$TOOLKIT_DIR/validate_bankit.py' \
                --tbl  '$OUT_DIR/bac7c17_gb/bac7c17_feature_table.tbl' \
                --nuc  '$OUT_DIR/bac7c17_gb/bac7c17_nucleotide.fasta' \
                --prot '$OUT_DIR/bac7c17_gb/bac7c17_protein.fasta' 2>&1"
    fi

    info "End-to-end output saved to: $OUT_DIR"
    echo ""
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo "========================================================================"
if [ "$FAIL_COUNT" -eq 0 ]; then
    success "All $PASS_COUNT tests passed!"
    echo "========================================================================"
    exit 0
else
    fail "$FAIL_COUNT test(s) FAILED. $PASS_COUNT test(s) passed."
    echo "========================================================================"
    exit 1
fi
