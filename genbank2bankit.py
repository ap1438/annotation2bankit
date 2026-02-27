#!/usr/bin/env python3
"""
genbank2bankit.py
=================
Extract gene annotations from existing GenBank, GFF, GFF3, or GTF files
and reformat them as NCBI BankIt-compatible submission files.

This is the REVERSE tool: it takes already-annotated files (e.g. from a
published GenBank record or an annotation pipeline) and re-generates the
three BankIt submission files.  Useful for:
  - Reformatting an existing GenBank record for resubmission
  - Extracting a subset of genes from a large annotation
  - Converting GFF3/GTF annotation to BankIt format
  - Preparing updated/corrected submissions

Output files
------------
1. <prefix>_feature_table.tbl   – 5-column BankIt feature table
2. <prefix>_nucleotide.fasta    – BankIt-formatted genomic FASTA
3. <prefix>_protein.fasta       – Protein sequences (for reference)

Supported input formats
-----------------------
  --format genbank   GenBank flat file (.gb, .genbank, .txt)
  --format gff3      GFF3 annotation (.gff3, .gff)
  --format gtf       GTF annotation (.gtf)
  --format gff2      GFF2 / generic GFF
  --format auto      Auto-detect (default)

Usage
-----
    # From GenBank flat file
    python genbank2bankit.py \\
        --input   BAC_7C17.gb \\
        --format  genbank \\
        --organism "Arabidopsis halleri subsp. halleri" \\
        --outdir  bankit_output/

    # From GFF3 + separate FASTA
    python genbank2bankit.py \\
        --input   BAC_7C17.gff3 \\
        --fasta   BAC_7C17_complete_fasta.fasta \\
        --organism "Arabidopsis halleri subsp. halleri" \\
        --outdir  bankit_output/

    # Extract only HMA4 genes
    python genbank2bankit.py \\
        --input     EU382073.gb \\
        --gene-list AhHMA4-1,AhHMA4-2 \\
        --organism  "Arabidopsis halleri subsp. halleri" \\
        --outdir    bankit_output/

Required arguments
------------------
  --input    FILE     Annotation file (GenBank / GFF3 / GTF / GFF2).
  --organism STR      Scientific organism name.

Optional arguments
------------------
  --fasta     FILE    Separate FASTA file (required for GFF if no embedded FASTA).
  --format    STR     Input format: auto|genbank|gff3|gtf|gff2 (default: auto).
  --outdir    DIR     Output directory (default: current directory).
  --prefix    STR     Output filename prefix (default: "output").
  --mol-type  STR     Molecule type (default: "genomic DNA").
  --isolate   STR     Isolate name.
  --strain    STR     Strain name.
  --clone     STR     Clone identifier.
  --country   STR     Country (NCBI format: "Country:Region").
  --date      STR     Collection date.
  --sub-species STR   Sub-species.
  --genotype  STR     Genotype.
  --gene-list STR     Comma-separated gene IDs or file path (one per line).
  --no-validate       Skip BankIt validation checks.
  --verbose           Print detailed progress.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib.genbank_parser   import parse_genbank
from lib.gff_parser       import parse_gff, detect_format
from lib.feature_writer   import (write_feature_table, write_nucleotide_fasta,
                                   write_protein_fasta)
from lib.fasta_utils      import read_fasta
from lib.models           import SequenceRecord, GeneModel
from lib.ncbi_validator   import run_all_validations

from typing import Dict, List, Optional, Set


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog        = 'genbank2bankit.py',
        description = 'Convert GenBank/GFF/GTF annotations to NCBI BankIt files.',
        formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    req = p.add_argument_group('Required')
    req.add_argument('--input',    required=True, metavar='FILE',
                     help='Annotation file (GenBank / GFF3 / GTF / GFF2).')
    req.add_argument('--organism', required=True, metavar='STR',
                     help='Scientific organism name.')

    opt = p.add_argument_group('Input options')
    opt.add_argument('--fasta',   default='', metavar='FILE',
                     help='Separate FASTA (needed for GFF without embedded FASTA).')
    opt.add_argument('--format',  default='auto',
                     choices=['auto', 'genbank', 'gff3', 'gtf', 'gff2'],
                     help='Input file format (default: auto-detect).')

    out = p.add_argument_group('Output')
    out.add_argument('--outdir',  default='.', metavar='DIR')
    out.add_argument('--prefix',  default='output', metavar='STR')

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

    flt = p.add_argument_group('Filtering')
    flt.add_argument('--gene-list', default='', metavar='STR_or_FILE',
                     help='Comma-separated gene names or file with one per line.')

    flg = p.add_argument_group('Flags')
    flg.add_argument('--no-validate', action='store_true')
    flg.add_argument('--verbose',     action='store_true')
    return p


# ─────────────────────────────────────────────────────────────────────────────
# Main logic
# ─────────────────────────────────────────────────────────────────────────────

def main(argv=None):
    parser = build_parser()
    args   = parser.parse_args(argv)

    sep = '─' * 70
    print(f'\n{sep}')
    print('genbank2bankit  –  GenBank/GFF/GTF → NCBI BankIt Converter')
    print(sep)

    if not os.path.isfile(args.input):
        sys.exit(f'ERROR: Input file not found: {args.input}')

    # ── Detect format ─────────────────────────────────────────────────────
    fmt = args.format
    if fmt == 'auto':
        fmt = _detect_format(args.input)
    print(f'\n[1/5] Input format detected: {fmt.upper()}')

    # ── Parse input ───────────────────────────────────────────────────────
    print(f'[2/5] Parsing annotation file: {args.input}')
    genomic_seqs: Dict[str, str] = {}
    gene_models:  List[GeneModel] = []
    seq_records:  List[SequenceRecord] = []

    if fmt == 'genbank':
        results = parse_genbank(args.input)
        for seq_rec, models in results:
            genomic_seqs[seq_rec.seq_id] = seq_rec.sequence
            # Inherit organism from GenBank if not overridden
            rec = SequenceRecord(
                seq_id      = seq_rec.seq_id,
                sequence    = seq_rec.sequence,
                organism    = args.organism or seq_rec.organism,
                mol_type    = args.mol_type  or seq_rec.mol_type,
                description = seq_rec.description,
                modifiers   = _merge_modifiers(seq_rec.modifiers,
                                               _build_modifiers(args)),
            )
            seq_records.append(rec)
            gene_models.extend(models)
        print(f'      Parsed {len(seq_records)} record(s), '
              f'{len(gene_models)} gene model(s).')

    else:
        # GFF3 / GTF / GFF2
        fasta_path = args.fasta if args.fasta else None
        seqs, models = parse_gff(args.input, fasta_path=fasta_path, fmt=fmt)
        genomic_seqs = seqs
        gene_models  = models
        print(f'      Parsed {len(gene_models)} gene model(s).')

    if not gene_models:
        sys.exit('ERROR: No gene models found in input file.')

    # ── Gene-list filter ──────────────────────────────────────────────────
    if args.gene_list:
        allowed = _parse_gene_list(args.gene_list)
        before  = len(gene_models)
        gene_models = [gm for gm in gene_models
                       if gm.gene_id in allowed
                       or gm.qualifiers.get('gene', '') in allowed
                       or gm.qualifiers.get('locus_tag', '') in allowed]
        print(f'[3/5] Gene-list filter: {before} → {len(gene_models)} model(s).')
    else:
        print(f'[3/5] No gene-list filter applied.')

    if not gene_models:
        sys.exit('ERROR: No gene models remain after filtering.')

    # ── Ensure product qualifier ──────────────────────────────────────────
    for gm in gene_models:
        if 'product' not in gm.qualifiers or not gm.qualifiers['product']:
            gm.qualifiers['product'] = 'hypothetical protein'
        if 'gene' not in gm.qualifiers and gm.gene_id:
            gm.qualifiers['gene'] = gm.gene_id

    # ── Build SequenceRecord list for GFF input ───────────────────────────
    if not seq_records:
        modifiers = _build_modifiers(args)
        needed_ids = list(dict.fromkeys(gm.seq_id for gm in gene_models))
        for sid in needed_ids:
            seq = genomic_seqs.get(sid, '')
            if not seq:
                print(f'  [WARN] No sequence found for "{sid}". '
                      f'Provide --fasta if needed.')
            rec = SequenceRecord(
                seq_id      = sid,
                sequence    = seq,
                organism    = args.organism,
                mol_type    = args.mol_type,
                description = 'genomic sequence',
                modifiers   = modifiers,
            )
            seq_records.append(rec)

    # ── Write output ──────────────────────────────────────────────────────
    os.makedirs(args.outdir, exist_ok=True)
    prefix   = os.path.join(args.outdir, args.prefix)
    tbl_path  = f'{prefix}_feature_table.tbl'
    nuc_path  = f'{prefix}_nucleotide.fasta'
    prot_path = f'{prefix}_protein.fasta'

    print(f'\n[4/5] Writing output files to: {args.outdir}/')
    write_feature_table(gene_models, tbl_path)
    write_nucleotide_fasta(seq_records, nuc_path)
    write_protein_fasta(gene_models, prot_path, genomic_seqs)

    print(f'      Feature table : {tbl_path}')
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
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _detect_format(filepath: str) -> str:
    """Detect input file format."""
    ext = os.path.splitext(filepath)[1].lower()
    if ext in ('.gb', '.genbank', '.gbk'):
        return 'genbank'
    if ext in ('.gff3',):
        return 'gff3'
    if ext in ('.gtf',):
        return 'gtf'
    if ext in ('.gff', '.gff2'):
        # Could be GFF2 or GFF3 – use content detection
        return detect_format(filepath)
    if ext in ('.txt',):
        # Check if it looks like GenBank flat file
        with open(filepath, 'r', errors='replace') as fh:
            first_line = fh.readline()
        if first_line.startswith('LOCUS'):
            return 'genbank'
    # Default to GFF3 content detection
    return detect_format(filepath)


def _parse_gene_list(arg: str) -> Set[str]:
    """Parse gene list from comma-separated string or file."""
    genes: Set[str] = set()
    if os.path.isfile(arg):
        with open(arg, 'r') as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith('#'):
                    genes.add(line)
    else:
        for g in arg.split(','):
            g = g.strip()
            if g:
                genes.add(g)
    return genes


def _build_modifiers(args) -> Dict[str, str]:
    modifiers = {}
    for key in ('isolate', 'strain', 'clone', 'country',
                 'collection_date', 'sub_species', 'genotype'):
        val = getattr(args, key, '')
        if val:
            modifiers[key] = val
    return modifiers


def _merge_modifiers(base: Dict[str, str],
                     override: Dict[str, str]) -> Dict[str, str]:
    """Merge two modifier dicts; override takes priority."""
    merged = dict(base)
    merged.update({k: v for k, v in override.items() if v})
    return merged


if __name__ == '__main__':
    sys.exit(main())
