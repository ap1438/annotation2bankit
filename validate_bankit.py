#!/usr/bin/env python3
"""
validate_bankit.py
==================
Standalone validation tool for NCBI BankIt submission files.

Validates:
  1. Feature table (.tbl)       – format, coordinates, required qualifiers
  2. Nucleotide FASTA           – format, BankIt modifiers, IUPAC chars
  3. Protein FASTA (optional)   – format, valid amino acids, start M
  4. Cross-file consistency     – SeqID matching, CDS length, translations

Usage
-----
    python validate_bankit.py \\
        --tbl    output_feature_table.tbl \\
        --nuc    output_nucleotide.fasta \\
        [--prot  output_protein.fasta] \\
        [--verbose] \\
        [--report output_validation_report.txt]

Arguments
---------
  --tbl    FILE  Feature table file (.tbl) – required.
  --nuc    FILE  Nucleotide FASTA file – required.
  --prot   FILE  Protein FASTA file – optional but recommended.
  --report FILE  Save validation report to this file.
  --verbose      Print detailed error messages.

Exit codes
----------
  0 – All checks passed (no errors). Files may be submitted to BankIt.
  1 – One or more errors found. Fix issues before submission.
  2 – A required file was not found.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib.ncbi_validator import (
    validate_feature_table,
    validate_nucleotide_fasta,
    validate_protein_fasta,
    validate_cross_consistency,
    run_all_validations,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog        = 'validate_bankit.py',
        description = 'Validate NCBI BankIt submission files.',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        epilog      = (
            'Examples:\n'
            '  python validate_bankit.py --tbl output.tbl --nuc output_nuc.fasta\n'
            '  python validate_bankit.py --tbl a.tbl --nuc a.fasta --prot a_prot.fasta '
            '--report validation.txt'
        ),
    )
    p.add_argument('--tbl',     required=True, metavar='FILE',
                   help='BankIt feature table (.tbl).')
    p.add_argument('--nuc',     required=True, metavar='FILE',
                   help='Nucleotide FASTA file.')
    p.add_argument('--prot',    default='',    metavar='FILE',
                   help='Protein FASTA file (optional).')
    p.add_argument('--report',  default='',    metavar='FILE',
                   help='Save report to file.')
    p.add_argument('--verbose', action='store_true',
                   help='Print detailed messages.')
    return p


def main(argv=None):
    parser = build_parser()
    args   = parser.parse_args(argv)

    # ── Check files exist ─────────────────────────────────────────────────
    for label, path in [('Feature table', args.tbl),
                         ('Nucleotide FASTA', args.nuc)]:
        if not os.path.isfile(path):
            print(f'ERROR: {label} not found: {path}', file=sys.stderr)
            return 2

    if args.prot and not os.path.isfile(args.prot):
        print(f'WARNING: Protein FASTA not found: {args.prot}', file=sys.stderr)
        args.prot = ''

    # ── Run validation ────────────────────────────────────────────────────
    passed = run_all_validations(
        tbl_path        = args.tbl,
        nuc_fasta_path  = args.nuc,
        prot_fasta_path = args.prot or None,
        verbose         = True,
    )

    # ── Save report ───────────────────────────────────────────────────────
    if args.report:
        import io
        import contextlib

        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            run_all_validations(
                tbl_path        = args.tbl,
                nuc_fasta_path  = args.nuc,
                prot_fasta_path = args.prot or None,
                verbose         = True,
            )
        report_text = buf.getvalue()
        with open(args.report, 'w') as fh:
            fh.write(report_text)
        print(f'\nValidation report saved to: {args.report}')

    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(main())
