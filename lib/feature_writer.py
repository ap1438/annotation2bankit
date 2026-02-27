"""
feature_writer.py
-----------------
Write NCBI BankIt-compatible output files:
  1. Feature table (.tbl)      – 5-column tab-delimited format
  2. Nucleotide FASTA          – with BankIt-formatted deflines
  3. Protein FASTA             – protein sequences for reference/tbl2asn

NCBI BankIt Feature Table Rules (implemented here)
---------------------------------------------------
Format reference:
  https://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html
  https://www.insdc.org/submitting-standards/feature-table-documentation/

1.  File begins with:  >Feature SeqID
    SeqID must exactly match the FASTA defline identifier (first token
    after '>').

2.  Feature line (5 columns, tab-delimited):
        start <TAB> stop <TAB> feature_key

3.  Continuation interval (multi-exon features):
        start <TAB> stop
    (no feature_key on continuation lines)

4.  Qualifier line:
        <TAB><TAB><TAB>qualifier_key <TAB> qualifier_value
    (3 leading tabs – empty start, stop, and feature_key columns)
    For boolean qualifiers (no value):
        <TAB><TAB><TAB>qualifier_key

5.  Strand convention:
      Plus  strand  – start < stop  (ascending)
      Minus strand  – start > stop  (5' end of the mRNA listed first)
    Multi-exon minus-strand features are listed in transcript order
    (5'→3'), so coordinates descend from first to last interval.

6.  Partial features:
      5'-partial: '<' prepended to the 5'-most coordinate
        (for plus strand: '<start' on the first interval)
        (for minus strand: '<stop' on the last interval's stop position,
         but in table format this means 'start>' for that interval...
         actually: for minus strand 5' partial, the symbol is '<' before
         the FIRST coordinate written, which is the highest genomic pos)
      3'-partial: '>' appended to the 3'-most coordinate
        (for plus strand: 'stop>' → written as '>stop' before the value)
        (for minus strand: the last written coordinate, the smallest
         genomic position, gets '<' → but NCBI uses '<' for 5' and '>'
         for 3'; for minus strand 3' partial, last interval stop < )

    NCBI actual rule (from the Feature Table document):
      '<' before a coordinate means "extends beyond the beginning of the
           sequence or before the feature begins".
      '>' before a coordinate means "extends beyond the end of the
           sequence or the feature end is unknown".
    For a minus-strand gene read 5'→3':
      5'-partial (no start codon) → the 5' end is unknown →
        on the first interval written (highest genomic), the first
        coordinate (the higher one) gets '<' to mean "5' partial".
      3'-partial (no stop codon) → on the last interval (lowest genomic),
        the second coordinate gets '>'.

    Implementation note:
      We represent partial using a boolean pair (is_partial_5, is_partial_3)
      and apply the bracket notation at write time based on strand.

7.  Feature hierarchy for a protein-coding gene:
        gene        – gene span (one interval)
        mRNA        – exon intervals (includes UTR)
        CDS         – coding exon intervals (no UTR)
    All three share the same qualifiers (gene name / locus_tag / product).

8.  Required qualifiers for BankIt genomic submissions:
      CDS:   product, (codon_start for partial 5' CDS)
      gene:  gene or locus_tag
      mRNA:  product

9.  The 'translation' qualifier is NOT included in the feature table;
    NCBI derives it from the CDS coordinates and the submitted sequence.
"""

from typing import Dict, List, Optional, Tuple
from .models import GeneModel, SequenceRecord
from .fasta_utils import (format_bankit_defline, write_fasta,
                          extract_cds_from_genome, translate_cds)


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

_TAB = '\t'
_Q   = '\t\t\t'   # 3-tab prefix for qualifier lines


def _coord_str(coord: int, partial_symbol: Optional[str] = None) -> str:
    """Format a coordinate, optionally prepending '<' or '>'."""
    if partial_symbol:
        return f'{partial_symbol}{coord}'
    return str(coord)


def _write_feature_intervals(
    fh,
    intervals: list,        # list of (start, stop) in TRANSCRIPT ORDER
    feature_key: str,
    is_partial_5: bool,
    is_partial_3: bool,
    strand: str,
) -> None:
    """
    Write one feature block (feature_key + all intervals) to fh.

    intervals must already be sorted in 5'→3' transcript order
    (ascending for +, descending for -), with each item being a
    (start, stop) tuple where start <= stop in genomic coords.

    For minus-strand output, we swap each pair so first > second.
    """
    for idx, (gstart, gstop) in enumerate(intervals):
        is_first = (idx == 0)
        is_last  = (idx == len(intervals) - 1)

        if strand == '+':
            col1 = gstart
            col2 = gstop
            p5_here = is_partial_5 and is_first
            p3_here = is_partial_3 and is_last
            c1 = _coord_str(col1, '<' if p5_here else None)
            c2 = _coord_str(col2, '>' if p3_here else None)
        else:
            # Minus strand: swap so higher coord is first (5' of mRNA)
            col1 = gstop
            col2 = gstart
            # 5'-partial: first written coordinate (highest) gets '<'
            p5_here = is_partial_5 and is_first
            # 3'-partial: last written coordinate (lowest) gets '>'
            p3_here = is_partial_3 and is_last
            c1 = _coord_str(col1, '<' if p5_here else None)
            c2 = _coord_str(col2, '>' if p3_here else None)

        if is_first:
            fh.write(f'{c1}{_TAB}{c2}{_TAB}{feature_key}\n')
        else:
            fh.write(f'{c1}{_TAB}{c2}\n')


def _write_qualifier(fh, key: str, value: Optional[str] = None) -> None:
    """Write one qualifier line (3 leading tabs)."""
    if value is not None:
        fh.write(f'{_Q}{key}{_TAB}{value}\n')
    else:
        fh.write(f'{_Q}{key}\n')


# ─────────────────────────────────────────────────────────────────────────────
# Feature table writer
# ─────────────────────────────────────────────────────────────────────────────

def write_feature_table(
    gene_models: List[GeneModel],
    output_path: str,
    include_misc_features: bool = False,
) -> None:
    """
    Write a BankIt-compatible feature table file.

    The file groups gene models by seq_id.  Within each seq_id block the
    features are sorted by their 5'-most genomic position.

    Parameters
    ----------
    gene_models           : list of GeneModel objects
    output_path           : path to write the .tbl file
    include_misc_features : if True, write misc_feature annotations stored
                            in the model's qualifiers under 'misc_feature'
                            key (not standard; use for debugging only)
    """
    # Group by seq_id, preserving insertion order
    by_seq: Dict[str, List[GeneModel]] = {}
    for gm in gene_models:
        by_seq.setdefault(gm.seq_id, []).append(gm)

    with open(output_path, 'w') as fh:
        for seq_id, models in by_seq.items():
            # Sort by 5'-most genomic coordinate
            models_sorted = sorted(
                models,
                key=lambda m: min(m.gene_start, m.gene_stop)
            )

            fh.write(f'>Feature {seq_id}\n')

            for gm in models_sorted:
                _write_gene_model(fh, gm)

    print(f'[feature_writer] Feature table written: {output_path}')


def _write_gene_model(fh, gm: GeneModel) -> None:
    """Write the gene, mRNA, and CDS features for one GeneModel."""

    strand = gm.strand

    # ── 1. gene feature ──────────────────────────────────────────────────
    if strand == '+':
        gene_c1 = _coord_str(gm.gene_start,
                              '<' if gm.is_partial_5 else None)
        gene_c2 = _coord_str(gm.gene_stop,
                              '>' if gm.is_partial_3 else None)
    else:
        gene_c1 = _coord_str(gm.gene_stop,
                              '<' if gm.is_partial_5 else None)
        gene_c2 = _coord_str(gm.gene_start,
                              '>' if gm.is_partial_3 else None)

    fh.write(f'{gene_c1}{_TAB}{gene_c2}{_TAB}gene\n')

    # gene/locus_tag qualifier
    if 'gene' in gm.qualifiers:
        _write_qualifier(fh, 'gene', gm.qualifiers['gene'])
    if 'locus_tag' in gm.qualifiers:
        _write_qualifier(fh, 'locus_tag', gm.qualifiers['locus_tag'])
    if 'pseudo' in gm.qualifiers:
        _write_qualifier(fh, 'pseudo')   # boolean qualifier

    # ── 2. mRNA feature ──────────────────────────────────────────────────
    # Use explicit exon intervals if available; fall back to CDS intervals
    # (acceptable for AUGUSTUS predictions without UTR annotation)
    if gm.exon_intervals:
        mrna_ivs_raw = [(iv.start, iv.stop)
                        for iv in gm.exon_sorted_for_table()]
    else:
        mrna_ivs_raw = [(iv.start, iv.stop)
                        for iv in gm.cds_sorted_for_table()]

    _write_feature_intervals(
        fh, mrna_ivs_raw, 'mRNA',
        gm.is_partial_5, gm.is_partial_3, strand,
    )
    if 'product' in gm.qualifiers:
        _write_qualifier(fh, 'product', gm.qualifiers['product'])
    if 'gene' in gm.qualifiers:
        _write_qualifier(fh, 'gene', gm.qualifiers['gene'])
    if 'locus_tag' in gm.qualifiers:
        _write_qualifier(fh, 'locus_tag', gm.qualifiers['locus_tag'])

    # ── 3. CDS feature ───────────────────────────────────────────────────
    cds_ivs_raw = [(iv.start, iv.stop)
                   for iv in gm.cds_sorted_for_table()]

    _write_feature_intervals(
        fh, cds_ivs_raw, 'CDS',
        gm.is_partial_5, gm.is_partial_3, strand,
    )

    if 'product' in gm.qualifiers:
        _write_qualifier(fh, 'product', gm.qualifiers['product'])
    if 'gene' in gm.qualifiers:
        _write_qualifier(fh, 'gene', gm.qualifiers['gene'])
    if 'locus_tag' in gm.qualifiers:
        _write_qualifier(fh, 'locus_tag', gm.qualifiers['locus_tag'])
    if 'protein_id' in gm.qualifiers:
        _write_qualifier(fh, 'protein_id', gm.qualifiers['protein_id'])
    if 'codon_start' in gm.qualifiers:
        _write_qualifier(fh, 'codon_start', gm.qualifiers['codon_start'])
    if 'note' in gm.qualifiers:
        _write_qualifier(fh, 'note', gm.qualifiers['note'])
    if 'pseudo' in gm.qualifiers:
        _write_qualifier(fh, 'pseudo')

    # ── Blank line between genes (improves readability) ──────────────────
    fh.write('\n')


# ─────────────────────────────────────────────────────────────────────────────
# Nucleotide FASTA writer
# ─────────────────────────────────────────────────────────────────────────────

def write_nucleotide_fasta(
    records: List[SequenceRecord],
    output_path: str,
    width: int = 60,
) -> None:
    """
    Write BankIt-formatted nucleotide FASTA file.

    Each SequenceRecord provides the seq_id, sequence, organism, mol_type,
    and optional source modifier key-value pairs.

    The defline format follows NCBI BankIt rules:
        >SeqID [organism=...] [mol_type=...] [modifier=value ...] description

    Parameters
    ----------
    records     : list of SequenceRecord
    output_path : path to write
    width       : nucleotide characters per line (default 60)
    """
    with open(output_path, 'w') as fh:
        for rec in records:
            defline = format_bankit_defline(
                seq_id      = rec.seq_id,
                organism    = rec.organism,
                mol_type    = rec.mol_type,
                description = rec.description,
                **rec.modifiers,
            )
            # defline already contains the '>' prefix
            fh.write(defline + '\n')
            seq = rec.sequence.upper()
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + '\n')

    print(f'[feature_writer] Nucleotide FASTA written: {output_path}')


# ─────────────────────────────────────────────────────────────────────────────
# Protein FASTA writer
# ─────────────────────────────────────────────────────────────────────────────

def write_protein_fasta(
    gene_models: List[GeneModel],
    output_path: str,
    genomic_seqs: Optional[Dict[str, str]] = None,
    width: int = 60,
) -> None:
    """
    Write a protein FASTA file for all gene models.

    Protein sequences are sourced (in priority order) from:
      1. gm.protein_seq  (set by AUGUSTUS parser from HTML inline sequences)
      2. Translated from gm.coding_seq (if available)
      3. Extracted from genomic_seqs and translated (if provided)

    The defline format:
        >SeqID_gene_id [gene=GeneName] [protein=product description]

    Parameters
    ----------
    gene_models  : list of GeneModel
    output_path  : path to write
    genomic_seqs : dict {seq_id: full_sequence} – used as fallback source
    width        : amino acids per line (default 60)
    """
    with open(output_path, 'w') as fh:
        written = 0
        for gm in gene_models:
            protein = _get_protein(gm, genomic_seqs)
            if protein is None:
                print(f'  [WARN] No protein sequence available for '
                      f'{gm.gene_id} ({gm.transcript_id}). Skipping.')
                continue

            # Build defline
            prot_id = f'{gm.seq_id}_{gm.gene_id}'
            parts   = [f'>{prot_id}']
            gene_name = gm.qualifiers.get('gene') or gm.gene_id
            if gene_name:
                parts.append(f'[gene={gene_name}]')
            product = gm.qualifiers.get('product', '')
            if product:
                parts.append(f'[protein={product}]')
            parts.append(f'[transcript_id={gm.transcript_id}]')
            fh.write(' '.join(parts) + '\n')

            protein_clean = protein.rstrip('*')
            for i in range(0, len(protein_clean), width):
                fh.write(protein_clean[i:i + width] + '\n')
            written += 1

    print(f'[feature_writer] Protein FASTA written: {output_path} '
          f'({written} sequences)')


def _get_protein(
    gm: GeneModel,
    genomic_seqs: Optional[Dict[str, str]],
) -> Optional[str]:
    """Return protein sequence for a GeneModel (see write_protein_fasta)."""
    if gm.protein_seq:
        return gm.protein_seq.rstrip('*')

    if gm.coding_seq:
        codon_start = int(gm.qualifiers.get('codon_start', 1))
        return translate_cds(gm.coding_seq, codon_start)

    if genomic_seqs and gm.seq_id in genomic_seqs:
        ivs = [(iv.start, iv.stop) for iv in gm.cds_intervals]
        try:
            cds_seq = extract_cds_from_genome(
                genomic_seqs[gm.seq_id], ivs, gm.strand)
            codon_start = int(gm.qualifiers.get('codon_start', 1))
            return translate_cds(cds_seq, codon_start)
        except Exception as e:
            print(f'  [WARN] Could not extract/translate CDS for '
                  f'{gm.gene_id}: {e}')

    return None
