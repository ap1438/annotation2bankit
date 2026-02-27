#!/usr/bin/env python3
"""
augustus2bankit.py
==================
Convert AUGUSTUS gene prediction output (HTML or plain GFF text) into
NCBI BankIt-compatible submission files.

Output files
------------
1. <prefix>_feature_table.tbl   – 5-column BankIt feature table
2. <prefix>_nucleotide.fasta    – BankIt-formatted genomic FASTA
3. <prefix>_protein.fasta       – Protein sequences (for reference)

Usage
-----
    python augustus2bankit.py \\
        --augustus  augustus_output.html \\
        --fasta     genomic_sequences.fasta \\
        --outdir    output/ \\
        --organism  "Arabidopsis halleri subsp. halleri" \\
        [OPTIONS]

Required arguments
------------------
  --augustus FILE     AUGUSTUS HTML or plain GFF output file.
  --fasta    FILE     Genomic FASTA file (source sequences submitted to AUGUSTUS).
  --organism STR      Scientific organism name (e.g. "Arabidopsis halleri subsp. halleri").

Optional arguments
------------------
  --outdir   DIR      Output directory (default: current directory).
  --prefix   STR      Output file prefix (default: "output").
  --mol-type STR      Molecule type (default: "genomic DNA").
  --isolate  STR      Isolate / accession name.
  --strain   STR      Strain name.
  --clone    STR      Clone identifier.
  --country  STR      Country of origin (NCBI format: "Country:Region").
  --date     STR      Collection date (e.g. "2007" or "Jan-2007").
  --sub-species STR   Sub-species name.
  --genotype STR      Genotype.
  --gene-name STR     Gene name qualifier to apply to all models
                      (overrides AUGUSTUS gene IDs; useful for single-gene runs).
  --product  STR      Product description to apply to all CDS features.
  --locus-tag STR     Locus tag prefix (e.g. "AhHMA4"; appended with -N).
  --gene-list FILE    Plain-text file listing gene IDs (one per line) to
                      include; all other predictions are skipped.
  --primary-only      Keep only the highest-scoring isoform per gene.
  --no-validate       Skip post-generation BankIt validation checks.
  --verbose           Print detailed progress.

Examples
--------
# Basic usage – single-gene BAC submission
python augustus2bankit.py \\
    --augustus  Augustus_Wall10123.html \\
    --fasta     Wall10.HMA4-123.fasta \\
    --organism  "Arabidopsis halleri subsp. halleri" \\
    --isolate   "Wall_10_v1.0.0" \\
    --product   "Zn/Cd P(IB)-type ATPase" \\
    --outdir    bankit_output/

# With gene naming and locus tags
python augustus2bankit.py \\
    --augustus  Augustus_result_display.html \\
    --fasta     BAC_7C17_complete_fasta.fasta \\
    --organism  "Arabidopsis halleri subsp. halleri" \\
    --clone     "BAC_7C17" \\
    --outdir    bankit_output/ \\
    --primary-only

NCBI BankIt submission notes
-----------------------------
After generating the files, upload them at:
  https://www.ncbi.nlm.nih.gov/WebSub/

Upload order:
  1. Nucleotide FASTA (one sequence per submitted entry)
  2. Feature table (.tbl) paired with FASTA
  3. Review annotations and confirm submission
"""

import argparse
import os
import sys

# ── Resolve lib/ directory regardless of working directory ────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib.augustus_parser import parse_augustus_output, filter_primary_isoforms
from lib.feature_writer   import (write_feature_table, write_nucleotide_fasta,
                                   write_protein_fasta)
from lib.fasta_utils      import read_fasta
from lib.models           import SequenceRecord, GeneModel
from lib.ncbi_validator   import run_all_validations

from typing import Dict, List, Optional


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog        = 'augustus2bankit.py',
        description = 'Convert AUGUSTUS output to NCBI BankIt submission files.',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        epilog      = __doc__.split('NCBI BankIt submission notes')[1],
    )
    req = p.add_argument_group('Required')
    req.add_argument('--augustus',  required=True, metavar='FILE',
                     help='AUGUSTUS HTML or plain GFF output file.')
    req.add_argument('--fasta',     required=True, metavar='FILE',
                     help='Genomic FASTA file (input to AUGUSTUS).')
    req.add_argument('--organism',  required=True, metavar='STR',
                     help='Scientific organism name.')

    opt = p.add_argument_group('Output')
    opt.add_argument('--outdir',    default='.', metavar='DIR',
                     help='Output directory (default: current dir).')
    opt.add_argument('--prefix',    default='output', metavar='STR',
                     help='Output filename prefix (default: output).')

    src = p.add_argument_group('Source modifiers (BankIt)')
    src.add_argument('--mol-type',    default='genomic DNA')
    src.add_argument('--isolate',     default='', metavar='STR')
    src.add_argument('--strain',      default='', metavar='STR')
    src.add_argument('--clone',       default='', metavar='STR')
    src.add_argument('--country',     default='', metavar='STR')
    src.add_argument('--date',        default='', metavar='STR',
                     dest='collection_date')
    src.add_argument('--sub-species', default='', metavar='STR',
                     dest='sub_species')
    src.add_argument('--genotype',    default='', metavar='STR')

    ann = p.add_argument_group('Annotation overrides')
    ann.add_argument('--gene-name',  default='', metavar='STR',
                     help='Override gene name for all models.')
    ann.add_argument('--product',    default='', metavar='STR',
                     help='Override product qualifier for all CDS.')
    ann.add_argument('--locus-tag',  default='', metavar='STR',
                     help='Locus tag prefix (e.g. AhHMA4); models get -1, -2 ...')
    ann.add_argument('--gene-list',  default='', metavar='FILE',
                     help='File with gene IDs to include (one per line).')

    flg = p.add_argument_group('Flags')
    flg.add_argument('--primary-only', action='store_true',
                     help='Keep only first isoform per gene.')
    flg.add_argument('--no-validate',  action='store_true',
                     help='Skip BankIt validation checks.')
    flg.add_argument('--verbose',      action='store_true')
    return p


# ─────────────────────────────────────────────────────────────────────────────
# Main logic
# ─────────────────────────────────────────────────────────────────────────────

def main(argv=None):
    parser = build_parser()
    args   = parser.parse_args(argv)

    sep = '─' * 70
    print(f'\n{sep}')
    print('augustus2bankit  –  AUGUSTUS → NCBI BankIt Converter')
    print(sep)

    # ── Load genomic FASTA ────────────────────────────────────────────────
    print(f'\n[1/5] Loading genomic FASTA: {args.fasta}')
    if not os.path.isfile(args.fasta):
        sys.exit(f'ERROR: FASTA file not found: {args.fasta}')
    genomic_seqs = read_fasta(args.fasta)
    print(f'      Loaded {len(genomic_seqs)} sequence(s): '
          f'{", ".join(list(genomic_seqs.keys())[:5])}')

    # ── Parse AUGUSTUS output ─────────────────────────────────────────────
    print(f'\n[2/5] Parsing AUGUSTUS output: {args.augustus}')
    if not os.path.isfile(args.augustus):
        sys.exit(f'ERROR: AUGUSTUS file not found: {args.augustus}')
    gene_models = parse_augustus_output(args.augustus)
    print(f'      Parsed {len(gene_models)} transcript model(s).')

    if not gene_models:
        sys.exit('ERROR: No gene models found in AUGUSTUS output. '
                 'Check that the file contains GFF predictions.')

    # ── Filter ───────────────────────────────────────────────────────────
    if args.primary_only:
        gene_models = filter_primary_isoforms(gene_models)
        print(f'      After --primary-only filter: {len(gene_models)} model(s).')

    if args.gene_list:
        allowed = _load_gene_list(args.gene_list)
        gene_models = [gm for gm in gene_models if gm.gene_id in allowed]
        print(f'      After --gene-list filter: {len(gene_models)} model(s).')

    if not gene_models:
        sys.exit('ERROR: No gene models remain after filtering.')

    # ── Apply annotation overrides ────────────────────────────────────────
    print(f'\n[3/5] Applying annotation overrides and qualifiers...')
    gene_models = _apply_overrides(gene_models, args)

    # ── Build SequenceRecord objects for nucleotide FASTA ─────────────────
    # Each unique seq_id referenced by gene models needs a record.
    modifiers = _build_modifiers(args)
    seq_ids_needed = list(dict.fromkeys(gm.seq_id for gm in gene_models))

    records: List[SequenceRecord] = []
    for sid in seq_ids_needed:
        if sid not in genomic_seqs:
            print(f'  [WARN] Sequence "{sid}" referenced in AUGUSTUS output '
                  f'but not found in FASTA. Skipping.')
            continue
        rec = SequenceRecord(
            seq_id      = sid,
            sequence    = genomic_seqs[sid],
            organism    = args.organism,
            mol_type    = args.mol_type,
            description = 'genomic sequence',
            modifiers   = modifiers,
        )
        records.append(rec)

    # ── Write output files ────────────────────────────────────────────────
    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, args.prefix)

    tbl_path  = f'{prefix}_feature_table.tbl'
    nuc_path  = f'{prefix}_nucleotide.fasta'
    prot_path = f'{prefix}_protein.fasta'

    print(f'\n[4/5] Writing output files to: {args.outdir}/')

    write_feature_table(gene_models, tbl_path)
    write_nucleotide_fasta(records, nuc_path)
    write_protein_fasta(gene_models, prot_path, genomic_seqs)

    print(f'\n      Feature table : {tbl_path}')
    print(f'      Nucleotide    : {nuc_path}')
    print(f'      Protein       : {prot_path}')

    # ── Validate ──────────────────────────────────────────────────────────
    if not args.no_validate:
        print(f'\n[5/5] Running NCBI BankIt validation checks...')
        passed = run_all_validations(tbl_path, nuc_path, prot_path,
                                     verbose=True)
        if not passed:
            print('\nPlease fix errors before submitting to NCBI BankIt.')
            return 1
    else:
        print('\n[5/5] Validation skipped (--no-validate).')

    print(f'\n{sep}')
    print('Conversion complete. Files are ready for BankIt submission.')
    print(f'{sep}\n')
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Helper functions
# ─────────────────────────────────────────────────────────────────────────────

def _load_gene_list(filepath: str) -> set:
    gene_set = set()
    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith('#'):
                gene_set.add(line)
    return gene_set


def _build_modifiers(args) -> Dict[str, str]:
    modifiers = {}
    for key in ('isolate', 'strain', 'clone', 'country',
                 'collection_date', 'sub_species', 'genotype'):
        val = getattr(args, key, '')
        if val:
            modifiers[key] = val
    return modifiers


def _apply_overrides(
    gene_models: List[GeneModel],
    args,
) -> List[GeneModel]:
    """
    Apply user-provided annotation overrides to all gene models.

    Overrides applied (when provided):
      --gene-name : set qualifiers['gene'] for all models
      --product   : set qualifiers['product'] for all models
      --locus-tag : set qualifiers['locus_tag'] = prefix-N
    """
    for idx, gm in enumerate(gene_models, 1):
        if args.gene_name:
            gm.qualifiers['gene'] = args.gene_name
        elif 'gene' not in gm.qualifiers:
            # Use AUGUSTUS gene_id as gene name
            gm.qualifiers['gene'] = gm.gene_id

        if args.product:
            gm.qualifiers['product'] = args.product
        elif 'product' not in gm.qualifiers:
            gm.qualifiers['product'] = 'hypothetical protein'

        if args.locus_tag:
            gm.qualifiers['locus_tag'] = f'{args.locus_tag}-{idx}'

    return gene_models


if __name__ == '__main__':
    sys.exit(main())
