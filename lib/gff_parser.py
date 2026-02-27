"""
gff_parser.py
-------------
Parse GFF2, GFF3, and GTF annotation files and return GeneModel objects.

Format detection
----------------
The format is detected automatically by inspecting the file header and
column 9 attribute syntax:
  GFF3 : header line "##gff-version 3", attributes use "key=value;" syntax
  GTF  : attributes use 'key "value";' syntax (similar to GFF2)
  GFF2 : fallback (attributes may be free-text or key "value";)

Location convention
-------------------
All three formats use 1-based, closed [start, stop] intervals with
start <= stop.  The strand column ('+', '-', or '.') indicates direction.

GFF3 parent-child relationship
--------------------------------
  gene  → mRNA  → CDS / exon
Gene models are assembled by traversing the Parent= attribute chain.

GTF relationship
-----------------
  gene_id and transcript_id attributes link features.
  "gene" and "transcript" features may or may not be present.

Embedded FASTA
--------------
GFF3 files may contain an embedded FASTA section beginning with "##FASTA".
If present, the sequences are extracted automatically.
"""

import re
from typing import Dict, List, Optional, Tuple

from .models import GeneModel, CDSInterval, ExonInterval, SequenceRecord
from .fasta_utils import read_fasta


# ─────────────────────────────────────────────────────────────────────────────
# Format detection
# ─────────────────────────────────────────────────────────────────────────────

def detect_format(filepath: str) -> str:
    """
    Detect GFF variant.

    Returns one of: 'gff3', 'gtf', 'gff2'
    """
    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('##gff-version 3') or line.startswith('##gff-version3'):
                return 'gff3'
            if line.startswith('#'):
                continue
            # Check attribute column
            cols = line.split('\t')
            if len(cols) >= 9:
                attr = cols[8]
                if 'ID=' in attr or 'Parent=' in attr:
                    return 'gff3'
                if re.search(r'gene_id\s+"', attr) or re.search(r'transcript_id\s+"', attr):
                    return 'gtf'
            break
    return 'gff2'


# ─────────────────────────────────────────────────────────────────────────────
# Attribute parsers
# ─────────────────────────────────────────────────────────────────────────────

def _parse_gff3_attrs(attr_str: str) -> Dict[str, str]:
    """Parse GFF3 key=value; attribute string."""
    attrs: Dict[str, str] = {}
    for part in attr_str.split(';'):
        part = part.strip()
        if '=' in part:
            k, _, v = part.partition('=')
            # URL-decode common characters
            v = v.replace('%2C', ',').replace('%3B', ';').replace('%3D', '=')
            attrs[k.strip()] = v.strip()
    return attrs


def _parse_gtf_attrs(attr_str: str) -> Dict[str, str]:
    """Parse GTF/GFF2 'key "value";' attribute string."""
    attrs: Dict[str, str] = {}
    for m in re.finditer(r'(\w+)\s+"([^"]*)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs


# ─────────────────────────────────────────────────────────────────────────────
# Embedded FASTA extractor
# ─────────────────────────────────────────────────────────────────────────────

def _split_gff3_fasta(filepath: str) -> Tuple[List[str], Dict[str, str]]:
    """
    Split a GFF3 file into annotation lines and embedded FASTA sequences.

    Returns (annotation_lines, {seq_id: sequence}).
    """
    ann_lines: List[str] = []
    fasta_lines: List[str] = []
    in_fasta = False

    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.strip() == '##FASTA' or line.strip().startswith('##FASTA'):
                in_fasta = True
                continue
            if in_fasta:
                fasta_lines.append(line)
            else:
                ann_lines.append(line)

    sequences: Dict[str, str] = {}
    if fasta_lines:
        cur_id: Optional[str] = None
        cur_seq: List[str] = []
        for line in fasta_lines:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if cur_id:
                    sequences[cur_id] = ''.join(cur_seq).upper()
                cur_id = line[1:].split()[0]
                cur_seq = []
            elif cur_id:
                cur_seq.append(line.strip())
        if cur_id:
            sequences[cur_id] = ''.join(cur_seq).upper()

    return ann_lines, sequences


# ─────────────────────────────────────────────────────────────────────────────
# GFF3 parser
# ─────────────────────────────────────────────────────────────────────────────

def _parse_gff3(
    lines: List[str],
) -> List[GeneModel]:
    """
    Parse GFF3 annotation lines into GeneModel objects.

    Builds a parent→children map, then assembles gene models by
    walking gene → mRNA → CDS/exon hierarchies.
    """
    # Raw features: id → dict
    features_by_id:  Dict[str, dict]        = {}
    children_of:     Dict[str, List[str]]   = {}  # parent_id → [child_ids]
    gene_ids:        List[str]              = []

    for line in lines:
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        cols = line.split('\t')
        if len(cols) < 8:
            continue

        seq_id  = cols[0]
        feature = cols[2]
        try:
            start = int(cols[3])
            stop  = int(cols[4])
        except ValueError:
            continue
        strand    = cols[6] if cols[6] in ('+', '-') else '+'
        frame_str = cols[7]
        try:
            frame = int(frame_str) if frame_str not in ('.', '') else 0
        except ValueError:
            frame = 0

        attrs = _parse_gff3_attrs(cols[8] if len(cols) > 8 else '')
        feat_id  = attrs.get('ID', '')
        parent   = attrs.get('Parent', '')

        feat_rec = {
            'seq_id':   seq_id,
            'feature':  feature,
            'start':    start,
            'stop':     stop,
            'strand':   strand,
            'frame':    frame,
            'attrs':    attrs,
            'id':       feat_id,
            'parent':   parent,
        }

        if feat_id:
            features_by_id[feat_id] = feat_rec
        if parent:
            children_of.setdefault(parent, []).append(feat_id)
        if feature == 'gene':
            gene_ids.append(feat_id)

    # ── Assemble gene models ──────────────────────────────────────────────
    models: List[GeneModel] = []

    for gene_id in gene_ids:
        gene_rec = features_by_id.get(gene_id, {})
        if not gene_rec:
            continue

        gene_attrs   = gene_rec.get('attrs', {})
        gene_seq_id  = gene_rec['seq_id']
        gene_strand  = gene_rec['strand']
        gene_start   = gene_rec['start']
        gene_stop    = gene_rec['stop']

        gene_name = (gene_attrs.get('Name')
                     or gene_attrs.get('gene')
                     or gene_attrs.get('locus_tag')
                     or gene_id)

        # Iterate over mRNA children
        mrna_ids = children_of.get(gene_id, [])
        if not mrna_ids:
            # No mRNA children – check if CDS directly under gene
            mrna_ids = [gene_id]   # treat gene as pseudo-transcript

        for mrna_idx, mrna_id in enumerate(mrna_ids):
            mrna_rec   = features_by_id.get(mrna_id, gene_rec)
            mrna_attrs = mrna_rec.get('attrs', {})

            product  = mrna_attrs.get('product', gene_attrs.get('product', ''))
            locus_tag = (mrna_attrs.get('locus_tag')
                         or gene_attrs.get('locus_tag', ''))
            protein_id = mrna_attrs.get('protein_id', '')
            note       = mrna_attrs.get('note', gene_attrs.get('note', ''))
            is_pseudo  = ('pseudo' in mrna_attrs or 'pseudo' in gene_attrs
                          or mrna_attrs.get('pseudogene', '') != ''
                          or gene_attrs.get('pseudogene', '') != '')

            t_id = mrna_id

            # Collect CDS and exon children of this mRNA
            cds_intervals: List[CDSInterval]  = []
            exon_intervals: List[ExonInterval] = []

            child_ids = children_of.get(mrna_id, [])
            # Also try direct children of gene if mRNA == gene
            if mrna_id == gene_id:
                child_ids = children_of.get(gene_id, [])

            for child_id in child_ids:
                child = features_by_id.get(child_id, {})
                if not child:
                    continue
                ctype = child.get('feature', '')
                cs    = child.get('start', 0)
                ce    = child.get('stop', 0)
                cf    = child.get('frame', 0)
                if ctype == 'CDS':
                    cds_intervals.append(CDSInterval(cs, ce, cf))
                elif ctype == 'exon':
                    exon_intervals.append(ExonInterval(cs, ce))

            if not cds_intervals and not is_pseudo:
                continue  # skip non-coding features

            # Partiality: detect from GFF3 partial=true or start_range/end_range
            partial_attr = mrna_attrs.get('partial', gene_attrs.get('partial', 'false'))
            is_partial_5 = partial_attr == 'true' or 'start_range' in mrna_attrs
            is_partial_3 = partial_attr == 'true' or 'end_range' in mrna_attrs

            # For NCBI GFF3: start_range=.,N means 5' partial; end_range=N,. means 3' partial
            start_range = mrna_attrs.get('start_range', '')
            end_range   = mrna_attrs.get('end_range', '')
            if start_range.startswith('.'):
                is_partial_5 = True
            if end_range.endswith('.'):
                is_partial_3 = True

            qualifiers: Dict[str, str] = {}
            qualifiers['gene'] = gene_name
            if locus_tag:
                qualifiers['locus_tag'] = locus_tag
            if product:
                qualifiers['product'] = product
            if protein_id:
                qualifiers['protein_id'] = protein_id
            if note:
                qualifiers['note'] = note
            if is_pseudo:
                qualifiers['pseudo'] = ''

            model = GeneModel(
                gene_id        = gene_name,
                transcript_id  = t_id,
                seq_id         = gene_seq_id,
                strand         = gene_strand,
                gene_start     = gene_start,
                gene_stop      = gene_stop,
                cds_intervals  = cds_intervals,
                exon_intervals = exon_intervals,
                is_partial_5   = is_partial_5,
                is_partial_3   = is_partial_3,
                qualifiers     = qualifiers,
            )
            models.append(model)

    return models


# ─────────────────────────────────────────────────────────────────────────────
# GTF parser
# ─────────────────────────────────────────────────────────────────────────────

def _parse_gtf(lines: List[str]) -> List[GeneModel]:
    """
    Parse GTF annotation lines into GeneModel objects.

    GTF uses gene_id and transcript_id attributes to link features.
    """
    # transcript_id → accumulated data
    transcript_data: Dict[str, dict] = {}
    gene_spans: Dict[str, Tuple[str, int, int, str]] = {}  # gene_id → (seq_id, start, stop, strand)

    for line in lines:
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        cols = line.split('\t')
        if len(cols) < 8:
            continue

        seq_id  = cols[0]
        feature = cols[2]
        try:
            start = int(cols[3])
            stop  = int(cols[4])
        except ValueError:
            continue
        strand    = cols[6] if cols[6] in ('+', '-') else '+'
        frame_str = cols[7]
        try:
            frame = int(frame_str) if frame_str not in ('.', '') else 0
        except ValueError:
            frame = 0

        attrs    = _parse_gtf_attrs(cols[8] if len(cols) > 8 else '')
        gene_id  = attrs.get('gene_id', '')
        t_id     = attrs.get('transcript_id', gene_id)
        product  = attrs.get('product', attrs.get('transcript_name', ''))
        locus_tag = attrs.get('locus_tag', '')
        protein_id = attrs.get('protein_id', '')
        note = attrs.get('note', '')

        if feature == 'gene' and gene_id:
            gene_spans[gene_id] = (seq_id, start, stop, strand)
            continue

        if feature in ('transcript', 'mRNA') and t_id:
            if t_id not in transcript_data:
                transcript_data[t_id] = {
                    'seq_id': seq_id, 'strand': strand,
                    'gene_id': gene_id, 't_start': start, 't_stop': stop,
                    'cds': [], 'exons': [],
                    'product': product, 'locus_tag': locus_tag,
                    'protein_id': protein_id, 'note': note,
                }
            continue

        if t_id and t_id not in transcript_data:
            transcript_data[t_id] = {
                'seq_id': seq_id, 'strand': strand,
                'gene_id': gene_id, 't_start': start, 't_stop': stop,
                'cds': [], 'exons': [],
                'product': product, 'locus_tag': locus_tag,
                'protein_id': protein_id, 'note': note,
            }

        if not t_id:
            continue
        td = transcript_data[t_id]

        if feature == 'CDS':
            td['cds'].append(CDSInterval(start, stop, frame))
        elif feature == 'exon':
            td['exons'].append(ExonInterval(start, stop))

        # Update product/etc. from CDS line if not set
        if not td['product'] and product:
            td['product'] = product
        if not td['locus_tag'] and locus_tag:
            td['locus_tag'] = locus_tag

    # ── Build GeneModel objects ───────────────────────────────────────────
    models: List[GeneModel] = []
    for t_id, td in transcript_data.items():
        if not td['cds']:
            continue
        gene_id = td['gene_id']
        if gene_id in gene_spans:
            _, gs, ge, _ = gene_spans[gene_id]
        else:
            gs = min(iv.start for iv in td['cds'])
            ge = max(iv.stop  for iv in td['cds'])

        qualifiers: Dict[str, str] = {}
        qualifiers['gene'] = gene_id
        if td['locus_tag']:
            qualifiers['locus_tag'] = td['locus_tag']
        if td['product']:
            qualifiers['product'] = td['product']
        if td['protein_id']:
            qualifiers['protein_id'] = td['protein_id']
        if td['note']:
            qualifiers['note'] = td['note']

        model = GeneModel(
            gene_id        = gene_id,
            transcript_id  = t_id,
            seq_id         = td['seq_id'],
            strand         = td['strand'],
            gene_start     = gs,
            gene_stop      = ge,
            cds_intervals  = td['cds'],
            exon_intervals = td['exons'],
            qualifiers     = qualifiers,
        )
        models.append(model)

    return models


# ─────────────────────────────────────────────────────────────────────────────
# GFF2 / generic parser
# ─────────────────────────────────────────────────────────────────────────────

def _parse_gff2(lines: List[str]) -> List[GeneModel]:
    """
    Minimal GFF2 parser: group CDS features by gene name in attribute col.
    """
    # For GFF2, group CDS by seq_id+gene_name
    gene_cds: Dict[str, dict] = {}

    for line in lines:
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        cols = line.split('\t')
        if len(cols) < 8:
            continue

        seq_id  = cols[0]
        feature = cols[2]
        if feature not in ('CDS', 'exon', 'gene', 'mRNA'):
            continue
        try:
            start = int(cols[3])
            stop  = int(cols[4])
        except ValueError:
            continue
        strand = cols[6] if cols[6] in ('+', '-') else '+'
        frame_str = cols[7]
        try:
            frame = int(frame_str) if frame_str not in ('.', '') else 0
        except ValueError:
            frame = 0

        attr_str = cols[8] if len(cols) > 8 else ''
        # Try GTF-style first, then plain
        attrs = _parse_gtf_attrs(attr_str)
        gene_id = (attrs.get('gene_id')
                   or attrs.get('gene_name')
                   or re.sub(r'"', '', attr_str.split(';')[0]).strip())

        key = f'{seq_id}|{gene_id}'
        if key not in gene_cds:
            gene_cds[key] = {
                'seq_id': seq_id, 'strand': strand,
                'gene_id': gene_id, 'cds': [], 'exons': [],
                'gs': start, 'ge': stop,
            }
        gd = gene_cds[key]
        gd['gs'] = min(gd['gs'], start)
        gd['ge'] = max(gd['ge'], stop)

        if feature == 'CDS':
            gd['cds'].append(CDSInterval(start, stop, frame))
        elif feature == 'exon':
            gd['exons'].append(ExonInterval(start, stop))

    models: List[GeneModel] = []
    for key, gd in gene_cds.items():
        if not gd['cds']:
            continue
        qualifiers = {'gene': gd['gene_id']}
        model = GeneModel(
            gene_id        = gd['gene_id'],
            transcript_id  = gd['gene_id'] + '.t1',
            seq_id         = gd['seq_id'],
            strand         = gd['strand'],
            gene_start     = gd['gs'],
            gene_stop      = gd['ge'],
            cds_intervals  = gd['cds'],
            exon_intervals = gd['exons'],
            qualifiers     = qualifiers,
        )
        models.append(model)
    return models


# ─────────────────────────────────────────────────────────────────────────────
# Public entry point
# ─────────────────────────────────────────────────────────────────────────────

def parse_gff(
    filepath: str,
    fasta_path: Optional[str] = None,
    fmt: str = 'auto',
) -> Tuple[Dict[str, str], List[GeneModel]]:
    """
    Parse a GFF/GFF3/GTF file into gene models.

    Parameters
    ----------
    filepath   : Path to the GFF/GFF3/GTF annotation file.
    fasta_path : Optional path to a separate FASTA file.  If the GFF3 has
                 an embedded ##FASTA section, this is ignored.
    fmt        : Format hint: 'auto' (default), 'gff3', 'gtf', 'gff2'.

    Returns
    -------
    (sequences_dict, gene_models)
    sequences_dict : {seq_id: sequence}  – may be empty if no FASTA provided
    gene_models    : list of GeneModel
    """
    # ── Detect or use specified format ────────────────────────────────────
    if fmt == 'auto':
        fmt = detect_format(filepath)
    print(f'[gff_parser] Detected format: {fmt.upper()}')

    # ── Extract annotation lines and embedded FASTA (GFF3) ────────────────
    if fmt == 'gff3':
        ann_lines, embedded_seqs = _split_gff3_fasta(filepath)
    else:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
            ann_lines = fh.readlines()
        embedded_seqs: Dict[str, str] = {}

    # ── Load external FASTA if needed ─────────────────────────────────────
    sequences: Dict[str, str] = embedded_seqs
    if not sequences and fasta_path:
        sequences = read_fasta(fasta_path)
        print(f'[gff_parser] Loaded {len(sequences)} sequences from {fasta_path}')

    # ── Parse annotation ──────────────────────────────────────────────────
    if fmt == 'gff3':
        models = _parse_gff3(ann_lines)
    elif fmt == 'gtf':
        models = _parse_gtf(ann_lines)
    else:
        models = _parse_gff2(ann_lines)

    print(f'[gff_parser] Parsed {len(models)} gene model(s).')
    return sequences, models
