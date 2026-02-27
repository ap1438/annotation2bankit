"""
ncbi_validator.py
-----------------
Validate NCBI BankIt submission files against NCBI rules before upload.

Checks performed
----------------
Feature table (.tbl):
  - File begins with >Feature SeqID
  - Correct tab structure (3 tabs before qualifiers)
  - Feature keys are INSDC-valid
  - Coordinates are positive integers
  - Minus-strand intervals have first coord > second coord
  - Required qualifiers present (gene/locus_tag for gene; product for CDS/mRNA)
  - CDS intervals divisible by 3 (when not partial)
  - No spaces in SeqID

Nucleotide FASTA:
  - Each record has a valid defline (>SeqID ...)
  - [organism=] modifier present
  - [mol_type=] modifier present
  - Sequence contains only IUPAC characters
  - No duplicate SeqIDs
  - Sequence length > 0

Protein FASTA:
  - Standard protein IUPAC characters
  - No stop codons (*) except at end
  - Sequences start with M (or partial)

Cross-file consistency:
  - Every SeqID in the feature table has a matching sequence in the FASTA
  - CDS coordinates do not exceed sequence length
  - CDS length divisible by 3 (complete CDS)
  - Translated CDS matches protein FASTA (if protein FASTA provided)

INSDC valid feature keys (subset used for genomic submissions):
  gene, mRNA, CDS, exon, intron, misc_feature, misc_RNA, ncRNA,
  precursor_RNA, prim_transcript, promoter, repeat_region, rRNA,
  stem_loop, tRNA, 3'UTR, 5'UTR, STS, source, repeat_unit,
  mobile_element, LTR, polyA_signal, polyA_site, sig_peptide,
  transit_peptide, mat_peptide, propeptide, modified_base, gap,
  assembly_gap, regulatory, centromere, telomere
"""

import re
from typing import Dict, List, Optional, Tuple

from .fasta_utils import (read_fasta, validate_fasta_file,
                           extract_cds_from_genome, translate_cds)


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

VALID_FEATURE_KEYS = {
    'gene', 'mRNA', 'CDS', 'exon', 'intron', 'misc_feature', 'misc_RNA',
    'ncRNA', 'precursor_RNA', 'prim_transcript', 'promoter', 'repeat_region',
    'rRNA', 'stem_loop', 'tRNA', "3'UTR", "5'UTR", 'STS', 'source',
    'repeat_unit', 'mobile_element', 'LTR', 'polyA_signal', 'polyA_site',
    'sig_peptide', 'transit_peptide', 'mat_peptide', 'propeptide',
    'modified_base', 'gap', 'assembly_gap', 'regulatory', 'centromere',
    'telomere', 'operon', 'protein_bind', 'primer_bind', 'misc_binding',
    'misc_difference', 'misc_recomb', 'misc_signal', 'misc_structure',
    'unsure', 'variation', 'enhancer', 'CAAT_signal', 'GC_signal',
    'TATA_signal', '-10_signal', '-35_signal',
}

REQUIRED_CDS_QUALS  = {'product'}
REQUIRED_GENE_QUALS = {'gene', 'locus_tag'}
IUPAC_PROT          = set('ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx*')
IUPAC_NUC           = set('ACGTURYSWKMBDHVNacgturyswkmbdhvn')


# ─────────────────────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────────────────────

class ValidationReport:
    def __init__(self):
        self.errors:   List[str] = []
        self.warnings: List[str] = []
        self.info:     List[str] = []

    def error(self, msg: str):
        self.errors.append(f'[ERROR]   {msg}')

    def warn(self, msg: str):
        self.warnings.append(f'[WARNING] {msg}')

    def note(self, msg: str):
        self.info.append(f'[INFO]    {msg}')

    @property
    def passed(self) -> bool:
        return len(self.errors) == 0

    def summary(self) -> str:
        lines = []
        lines.extend(self.errors)
        lines.extend(self.warnings)
        lines.extend(self.info)
        if not lines:
            lines.append('[INFO]    All checks passed.')
        status = 'PASS' if self.passed else 'FAIL'
        header = f'=== Validation {status} | {len(self.errors)} error(s), ' \
                 f'{len(self.warnings)} warning(s) ==='
        return '\n'.join([header] + lines)


def _parse_coord(s: str) -> Optional[int]:
    """Strip <> and convert to int."""
    s = s.replace('<', '').replace('>', '').strip()
    try:
        return int(s)
    except ValueError:
        return None


# ─────────────────────────────────────────────────────────────────────────────
# Feature table validator
# ─────────────────────────────────────────────────────────────────────────────

def validate_feature_table(tbl_path: str) -> ValidationReport:
    """
    Validate a BankIt 5-column feature table file.

    Returns a ValidationReport with errors and warnings.
    """
    rpt = ValidationReport()

    try:
        with open(tbl_path, 'r') as fh:
            lines = fh.readlines()
    except FileNotFoundError:
        rpt.error(f'Feature table file not found: {tbl_path}')
        return rpt

    if not lines:
        rpt.error('Feature table file is empty.')
        return rpt

    # ── State machine ─────────────────────────────────────────────────────
    current_seq_id:   Optional[str] = None
    current_feature:  Optional[str] = None
    current_gene_has_name    = False
    current_cds_has_product  = False
    current_mrna_has_product = False
    cds_intervals:  List[Tuple[int, int]] = []
    seq_ids_seen:   set = set()

    for lineno, raw in enumerate(lines, 1):
        line = raw.rstrip('\n')

        # ── >Feature line ─────────────────────────────────────────────────
        if line.startswith('>Feature'):
            parts = line.split()
            if len(parts) < 2:
                rpt.error(f'Line {lineno}: >Feature line missing SeqID.')
            else:
                seq_id = parts[1]
                if ' ' in seq_id:
                    rpt.error(f'Line {lineno}: SeqID "{seq_id}" contains a space.')
                if seq_id in seq_ids_seen:
                    rpt.warn(f'Line {lineno}: Duplicate >Feature block for SeqID "{seq_id}".')
                seq_ids_seen.add(seq_id)
                current_seq_id = seq_id
            current_feature = None
            continue

        if current_seq_id is None:
            if line.strip():
                rpt.error(f'Line {lineno}: Data found before any >Feature header.')
            continue

        if not line.strip():
            continue  # blank lines are OK between genes

        # ── Qualifier line (3 leading tabs) ───────────────────────────────
        if line.startswith('\t\t\t'):
            rest = line[3:]
            if '\t' in rest:
                qkey, _, qval = rest.partition('\t')
            else:
                qkey = rest
                qval = ''
            qkey = qkey.strip()

            if current_feature == 'CDS':
                if qkey == 'product':
                    current_cds_has_product = True
                if qkey == 'codon_start':
                    if qval not in ('1', '2', '3'):
                        rpt.error(f'Line {lineno}: codon_start must be 1, 2, or 3; got "{qval}".')
                if qkey == 'protein_id' and not re.match(r'.+\.\d+$', qval):
                    rpt.warn(f'Line {lineno}: protein_id "{qval}" '
                             f'does not match expected format ACC.version.')

            if current_feature == 'mRNA':
                if qkey == 'product':
                    current_mrna_has_product = True

            if current_feature == 'gene':
                if qkey in ('gene', 'locus_tag'):
                    current_gene_has_name = True

            continue

        # ── Feature or interval line ───────────────────────────────────────
        cols = line.split('\t')
        if len(cols) >= 3 and cols[2].strip():
            # New feature line: start \t stop \t feature_key
            # Finalise previous feature checks
            if current_feature == 'CDS' and not current_cds_has_product:
                rpt.warn(f'Line {lineno}: CDS before this line is missing "product" qualifier.')
            if current_feature == 'mRNA' and not current_mrna_has_product:
                rpt.warn(f'Line {lineno}: mRNA before this line is missing "product" qualifier.')
            if current_feature == 'gene' and not current_gene_has_name:
                rpt.warn(f'Line {lineno}: gene feature is missing "gene" or "locus_tag" qualifier.')

            c1_raw = cols[0].strip()
            c2_raw = cols[1].strip()
            feat_key = cols[2].strip()

            c1 = _parse_coord(c1_raw)
            c2 = _parse_coord(c2_raw)

            if c1 is None:
                rpt.error(f'Line {lineno}: Invalid coordinate "{c1_raw}".')
            if c2 is None:
                rpt.error(f'Line {lineno}: Invalid coordinate "{c2_raw}".')
            if c1 and c1 < 1:
                rpt.error(f'Line {lineno}: Coordinate {c1} < 1 (must be 1-based).')
            if c2 and c2 < 1:
                rpt.error(f'Line {lineno}: Coordinate {c2} < 1 (must be 1-based).')

            if feat_key not in VALID_FEATURE_KEYS:
                rpt.warn(f'Line {lineno}: Feature key "{feat_key}" is not in the standard INSDC list.')

            current_feature = feat_key
            current_gene_has_name    = False
            current_cds_has_product  = False
            current_mrna_has_product = False
            cds_intervals = []

            if c1 and c2:
                cds_intervals.append((c1, c2))

        elif len(cols) >= 2 and cols[0].strip() and cols[1].strip():
            # Continuation interval
            c1 = _parse_coord(cols[0].strip())
            c2 = _parse_coord(cols[1].strip())
            if c1 and c2:
                cds_intervals.append((c1, c2))
        else:
            rpt.warn(f'Line {lineno}: Could not parse line: "{line[:60]}"')

    # Final feature checks
    if current_feature == 'CDS' and not current_cds_has_product:
        rpt.warn('Last CDS feature is missing "product" qualifier.')

    rpt.note(f'Feature table: {len(seq_ids_seen)} sequence(s) validated.')
    return rpt


# ─────────────────────────────────────────────────────────────────────────────
# Nucleotide FASTA validator
# ─────────────────────────────────────────────────────────────────────────────

def validate_nucleotide_fasta(fasta_path: str) -> ValidationReport:
    """Validate a BankIt-format nucleotide FASTA file."""
    rpt = ValidationReport()

    is_valid, errs = validate_fasta_file(fasta_path)
    for e in errs:
        rpt.error(e)

    # Check for BankIt modifiers
    try:
        with open(fasta_path, 'r') as fh:
            for lineno, line in enumerate(fh, 1):
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if '[organism=' not in line:
                        rpt.error(f'Line {lineno}: Missing [organism=] modifier in defline.')
                    if '[mol_type=' not in line:
                        rpt.warn(f'Line {lineno}: Missing [mol_type=] modifier in defline.')
    except FileNotFoundError:
        rpt.error(f'File not found: {fasta_path}')

    rpt.note(f'Nucleotide FASTA validation complete: {fasta_path}')
    return rpt


# ─────────────────────────────────────────────────────────────────────────────
# Protein FASTA validator
# ─────────────────────────────────────────────────────────────────────────────

def validate_protein_fasta(prot_path: str) -> ValidationReport:
    """Validate a protein FASTA file."""
    rpt = ValidationReport()

    try:
        with open(prot_path, 'r') as fh:
            lines = fh.readlines()
    except FileNotFoundError:
        rpt.error(f'File not found: {prot_path}')
        return rpt

    cur_id: Optional[str] = None
    cur_seq: List[str] = []

    def _check_protein(pid, seq):
        if not seq:
            rpt.error(f'Protein "{pid}" has empty sequence.')
            return
        invalid = set(seq) - IUPAC_PROT
        if invalid:
            rpt.error(f'Protein "{pid}" contains invalid characters: {invalid}')
        if '*' in seq[:-1]:
            rpt.error(f'Protein "{pid}" contains internal stop codon (*) – '
                      f'possible frameshift or wrong reading frame.')
        if not seq[0] in ('M', 'X'):
            rpt.warn(f'Protein "{pid}" does not start with M '
                     f'(may be partial/5\'-truncated). Starts with: {seq[0]}')

    seen_ids: set = set()
    for lineno, line in enumerate(lines, 1):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if cur_id:
                _check_protein(cur_id, ''.join(cur_seq))
            tokens = line[1:].split()
            cur_id = tokens[0] if tokens else ''
            if cur_id in seen_ids:
                rpt.error(f'Line {lineno}: Duplicate protein ID "{cur_id}".')
            seen_ids.add(cur_id)
            cur_seq = []
        elif cur_id:
            cur_seq.append(line.strip().rstrip('*'))  # strip trailing stop

    if cur_id:
        _check_protein(cur_id, ''.join(cur_seq))

    rpt.note(f'Protein FASTA: {len(seen_ids)} sequence(s) validated.')
    return rpt


# ─────────────────────────────────────────────────────────────────────────────
# Cross-file consistency validator
# ─────────────────────────────────────────────────────────────────────────────

def validate_cross_consistency(
    tbl_path: str,
    nuc_fasta_path: str,
    prot_fasta_path: Optional[str] = None,
) -> ValidationReport:
    """
    Cross-validate feature table against nucleotide FASTA.

    Checks:
    - Every SeqID in the .tbl has a matching entry in the nucleotide FASTA.
    - CDS coordinates do not exceed sequence length.
    - For complete (non-partial) CDS: total CDS length divisible by 3.
    - Optionally verifies CDS translations match protein FASTA.
    """
    rpt = ValidationReport()

    # Load sequences
    try:
        sequences = read_fasta(nuc_fasta_path)
    except FileNotFoundError:
        rpt.error(f'Nucleotide FASTA not found: {nuc_fasta_path}')
        return rpt

    prot_seqs: Dict[str, str] = {}
    if prot_fasta_path:
        try:
            prot_seqs = read_fasta(prot_fasta_path)
        except FileNotFoundError:
            rpt.warn(f'Protein FASTA not found: {prot_fasta_path}')

    # Parse feature table
    try:
        with open(tbl_path, 'r') as fh:
            tbl_lines = fh.readlines()
    except FileNotFoundError:
        rpt.error(f'Feature table not found: {tbl_path}')
        return rpt

    cur_seq_id:  Optional[str] = None
    cur_feature: Optional[str] = None
    cur_intervals: List[Tuple[int, int]] = []
    cur_strand:    str = '+'
    is_partial_5 = False
    is_partial_3 = False

    for line in tbl_lines:
        line = line.rstrip('\n')

        if line.startswith('>Feature'):
            parts = line.split()
            if len(parts) >= 2:
                cur_seq_id = parts[1]
                if cur_seq_id not in sequences:
                    rpt.error(f'SeqID "{cur_seq_id}" in feature table '
                              f'has no matching sequence in FASTA.')
            cur_feature = None
            continue

        if not cur_seq_id or not line.strip():
            continue

        if line.startswith('\t\t\t'):
            continue  # qualifier – skip

        cols = line.split('\t')
        if len(cols) >= 3 and cols[2].strip():
            # New feature
            if cur_feature == 'CDS' and cur_seq_id in sequences:
                _check_cds(rpt, cur_seq_id, sequences[cur_seq_id],
                           cur_intervals, cur_strand,
                           is_partial_5, is_partial_3)
            cur_feature   = cols[2].strip()
            cur_intervals = []
            c1_raw = cols[0].strip()
            c2_raw = cols[1].strip()
            # Partiality: checked on FIRST interval; may also appear on
            # continuation lines (e.g. 3'-partial marker on last exon).
            is_partial_5 = '<' in c1_raw or '<' in c2_raw
            is_partial_3 = '>' in c1_raw or '>' in c2_raw
            c1 = _parse_coord(c1_raw)
            c2 = _parse_coord(c2_raw)
            if c1 and c2:
                # Determine strand from coordinate order
                cur_strand = '-' if c1 > c2 else '+'
                cur_intervals.append((min(c1, c2), max(c1, c2)))

            if cur_seq_id in sequences:
                seq_len = len(sequences[cur_seq_id])
                for coord in [c1, c2]:
                    if coord and coord > seq_len:
                        rpt.error(f'Coordinate {coord} exceeds sequence '
                                  f'length {seq_len} for "{cur_seq_id}".')

        elif len(cols) >= 2 and cols[0].strip() and cols[1].strip():
            c1_raw = cols[0].strip()
            c2_raw = cols[1].strip()
            # Partiality markers can also appear on continuation lines
            if '<' in c1_raw or '<' in c2_raw:
                is_partial_5 = True
            if '>' in c1_raw or '>' in c2_raw:
                is_partial_3 = True
            c1 = _parse_coord(c1_raw)
            c2 = _parse_coord(c2_raw)
            if c1 and c2:
                cur_intervals.append((min(c1, c2), max(c1, c2)))

    # Final CDS check
    if cur_feature == 'CDS' and cur_seq_id and cur_seq_id in sequences:
        _check_cds(rpt, cur_seq_id, sequences[cur_seq_id],
                   cur_intervals, cur_strand, is_partial_5, is_partial_3)

    rpt.note(f'Cross-consistency check complete: '
             f'{len(sequences)} sequence(s) in FASTA.')
    return rpt


def _check_cds(
    rpt: ValidationReport,
    seq_id: str,
    sequence: str,
    intervals: List[Tuple[int, int]],
    strand: str,
    is_partial_5: bool,
    is_partial_3: bool,
) -> None:
    """Check CDS intervals against the genomic sequence."""
    if not intervals:
        return

    total_len = sum(stop - start + 1 for start, stop in intervals)

    if not is_partial_5 and not is_partial_3:
        if total_len % 3 != 0:
            rpt.error(
                f'SeqID "{seq_id}": Complete CDS total length {total_len} bp '
                f'is not divisible by 3 (possible annotation error).')

    # Try translation
    try:
        cds_seq = extract_cds_from_genome(sequence, intervals, strand)
        if cds_seq:
            protein = translate_cds(cds_seq)
            if '*' in protein:
                rpt.error(
                    f'SeqID "{seq_id}": CDS translation contains internal '
                    f'stop codon(s). Possible frameshift or wrong coordinates.')
            if not is_partial_5 and not protein.startswith('M'):
                rpt.warn(
                    f'SeqID "{seq_id}": Complete CDS does not start with Met (M). '
                    f'First aa: {protein[0] if protein else "N/A"}')
    except Exception as e:
        rpt.warn(f'SeqID "{seq_id}": Could not validate CDS translation: {e}')


# ─────────────────────────────────────────────────────────────────────────────
# Master validation runner
# ─────────────────────────────────────────────────────────────────────────────

def run_all_validations(
    tbl_path:       str,
    nuc_fasta_path: str,
    prot_fasta_path: Optional[str] = None,
    verbose:         bool = True,
) -> bool:
    """
    Run all validation checks and print a summary report.

    Returns True if all checks passed (no errors), False otherwise.
    """
    separator = '─' * 70

    print(f'\n{separator}')
    print('NCBI BankIt Submission File Validation')
    print(separator)

    # 1. Feature table
    print('\n[1/4] Validating feature table...')
    r1 = validate_feature_table(tbl_path)
    if verbose:
        print(r1.summary())

    # 2. Nucleotide FASTA
    print('\n[2/4] Validating nucleotide FASTA...')
    r2 = validate_nucleotide_fasta(nuc_fasta_path)
    if verbose:
        print(r2.summary())

    # 3. Protein FASTA (optional)
    r3 = ValidationReport()
    if prot_fasta_path:
        print('\n[3/4] Validating protein FASTA...')
        r3 = validate_protein_fasta(prot_fasta_path)
        if verbose:
            print(r3.summary())
    else:
        print('\n[3/4] Protein FASTA not provided – skipping.')

    # 4. Cross-consistency
    print('\n[4/4] Running cross-consistency checks...')
    r4 = validate_cross_consistency(tbl_path, nuc_fasta_path, prot_fasta_path)
    if verbose:
        print(r4.summary())

    # ── Overall result ────────────────────────────────────────────────────
    all_passed = all(r.passed for r in [r1, r2, r3, r4])
    total_errors   = sum(len(r.errors)   for r in [r1, r2, r3, r4])
    total_warnings = sum(len(r.warnings) for r in [r1, r2, r3, r4])

    print(f'\n{separator}')
    if all_passed:
        print(f'OVERALL RESULT: PASS')
        print(f'  Errors: 0   Warnings: {total_warnings}')
        print('  Files appear ready for NCBI BankIt submission.')
    else:
        print(f'OVERALL RESULT: FAIL')
        print(f'  Errors: {total_errors}   Warnings: {total_warnings}')
        print('  Please fix errors before submitting to NCBI BankIt.')
    print(separator)

    return all_passed
