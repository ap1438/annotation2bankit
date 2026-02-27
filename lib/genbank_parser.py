"""
genbank_parser.py
-----------------
Parse NCBI GenBank flat file format (.gb / .genbank / .txt).

Extracts gene models (gene, mRNA, CDS features) and sequences.

GenBank flat file structure (abbreviated):
  LOCUS       name   length bp  DNA  linear   division  date
  DEFINITION  description
  ACCESSION   accession
  VERSION     version
  FEATURES             Location/Qualifiers
       source   1..total_length
                /organism="..."
                /mol_type="genomic DNA"
       gene     complement(100..500)
                /gene="GeneName"
       mRNA     join(100..200,300..500)
                /product="protein name"
       CDS      join(100..200,300..500)
                /gene="GeneName"
                /product="protein name"
                /translation="MAAA..."
  ORIGIN
          1 acgtacgt...
  //

Location format (handled here):
  Simple         : 100..500
  Complement     : complement(100..500)
  Join           : join(100..200,300..500)
  Complement join: complement(join(100..200,300..500))
  Partial        : <100..500  (5' partial) or 100..>500 (3' partial)
  Single base    : 100
  Between        : 100^101
"""

import re
from typing import Dict, List, Optional, Tuple

from .models import GeneModel, CDSInterval, ExonInterval, SequenceRecord
from .fasta_utils import extract_cds_from_genome, translate_cds


# ─────────────────────────────────────────────────────────────────────────────
# Location string parser
# ─────────────────────────────────────────────────────────────────────────────

def _parse_location(loc_str: str) -> Tuple[List[Tuple[int, int]], str, bool, bool]:
    """
    Parse a GenBank location string.

    Returns
    -------
    intervals   : list of (start, stop) tuples – 1-based, start<=stop
    strand      : '+' or '-'
    is_partial_5: True if '<' found on 5'-most coordinate
    is_partial_3: True if '>' found on 3'-most coordinate
    """
    loc = loc_str.strip()
    strand = '+'

    if loc.startswith('complement(') and loc.endswith(')'):
        strand = '-'
        loc = loc[len('complement('):-1].strip()

    # Handle complement(join(...))
    if loc.startswith('join(') and loc.endswith(')'):
        loc = loc[len('join('):-1]

    # Detect partiality
    is_partial_5 = '<' in loc
    is_partial_3 = '>' in loc

    # Remove partial markers for coordinate parsing
    loc_clean = loc.replace('<', '').replace('>', '')

    intervals: List[Tuple[int, int]] = []
    for part in loc_clean.split(','):
        part = part.strip()
        if '..' in part:
            s, e = part.split('..')
            try:
                intervals.append((int(s.strip()), int(e.strip())))
            except ValueError:
                continue
        elif part.isdigit():
            pos = int(part)
            intervals.append((pos, pos))
        elif '^' in part:
            # Between two positions – skip (insertion sites)
            pass

    return intervals, strand, is_partial_5, is_partial_3


# ─────────────────────────────────────────────────────────────────────────────
# GenBank flat file parser
# ─────────────────────────────────────────────────────────────────────────────

class _FeatureEntry:
    """Temporary holder during GenBank parsing."""
    def __init__(self, key: str, location: str):
        self.key        = key
        self.location   = location
        self.qualifiers: Dict[str, str] = {}
        self._cur_qual: Optional[str]   = None
        self._cur_val: List[str]        = []

    def finalise_qualifier(self):
        if self._cur_qual:
            self.qualifiers[self._cur_qual] = (
                ''.join(self._cur_val)
                .replace('"', '')
                .strip()
            )
        self._cur_qual = None
        self._cur_val  = []


def parse_genbank(filepath: str) -> List[Tuple[SequenceRecord, List[GeneModel]]]:
    """
    Parse a GenBank flat file and return a list of
    (SequenceRecord, [GeneModel, ...]) tuples – one per record (//  block).

    Multiple records in one file (multi-record GenBank) are fully supported.
    """
    results: List[Tuple[SequenceRecord, List[GeneModel]]] = []

    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        record_lines = []
        # Collect lines until '//'
        while i < len(lines) and not lines[i].startswith('//'):
            record_lines.append(lines[i])
            i += 1
        i += 1  # skip '//'

        if not record_lines:
            continue

        rec, models = _parse_genbank_record(record_lines)
        if rec is not None:
            results.append((rec, models))

    return results


def _parse_genbank_record(
    lines: List[str],
) -> Tuple[Optional[SequenceRecord], List[GeneModel]]:
    """Parse one GenBank record (everything up to //)."""

    # ── Pass 1: collect metadata and features ────────────────────────────
    locus_name:  str = ''
    seq_length:  int = 0
    definition:  str = ''
    accession:   str = ''
    organism:    str = ''

    section = 'header'
    features: List[_FeatureEntry] = []
    cur_feat: Optional[_FeatureEntry] = None

    seq_lines: List[str] = []
    in_origin = False

    for line in lines:
        if in_origin:
            # ORIGIN section – collect sequence digits
            seq_lines.append(
                re.sub(r'[^acgtnACGTN]', '', line))
            continue

        if line.startswith('LOCUS'):
            parts = line.split()
            if len(parts) >= 3:
                locus_name = parts[1]
                try:
                    seq_length = int(parts[2])
                except ValueError:
                    pass
            section = 'header'
            continue

        if line.startswith('DEFINITION'):
            definition = line[12:].rstrip()
            section = 'definition'
            continue

        if line.startswith('ACCESSION'):
            accession = line[12:].split()[0] if len(line) > 12 else ''
            section = 'header'
            continue

        if line.startswith('  ORGANISM'):
            organism = line[12:].rstrip()
            section = 'header'
            continue

        if line.startswith('FEATURES'):
            if cur_feat:
                cur_feat.finalise_qualifier()
                features.append(cur_feat)
                cur_feat = None
            section = 'features'
            continue

        if line.startswith('ORIGIN'):
            if cur_feat:
                cur_feat.finalise_qualifier()
                features.append(cur_feat)
                cur_feat = None
            section = 'origin'
            in_origin = True
            continue

        if section == 'definition' and line.startswith(' ' * 12):
            definition += ' ' + line.strip()
            continue

        if section == 'features':
            _process_feature_line(line, features, cur_feat)
            # Rebuild cur_feat reference after potential mutation
            if features and line.startswith('     ') and not line.startswith('      ' * 2):
                feat_line_match = re.match(r'     (\S+)\s+(\S.*)$', line)
                if feat_line_match:
                    if cur_feat:
                        cur_feat.finalise_qualifier()
                    cur_feat = _FeatureEntry(
                        feat_line_match.group(1),
                        feat_line_match.group(2).strip(),
                    )
                    # Replace the one we just added via _process_feature_line
                    # (we process twice – simplify by doing it in one pass)
            continue

    # Final feature
    if cur_feat:
        cur_feat.finalise_qualifier()
        features.append(cur_feat)

    # Rebuild with single-pass approach
    features, source_qual = _parse_features_singlepass(
        [l for l in lines
         if not l.startswith('ORIGIN') and not l.startswith('//')]
    )

    # ── Reconstruct sequence ──────────────────────────────────────────────
    sequence = ''.join(seq_lines).upper()

    # ── Build SequenceRecord ──────────────────────────────────────────────
    seq_id = accession or locus_name or 'UNKNOWN'
    rec = SequenceRecord(
        seq_id      = seq_id,
        sequence    = sequence,
        organism    = source_qual.get('organism', organism or 'Unknown'),
        mol_type    = source_qual.get('mol_type', 'genomic DNA'),
        description = definition.strip() or 'genomic sequence',
        modifiers   = {
            k: v for k, v in source_qual.items()
            if k not in ('organism', 'mol_type')
        },
    )

    # ── Build GeneModel objects ───────────────────────────────────────────
    models = _build_gene_models(features, seq_id)

    return rec, models


def _parse_features_singlepass(
    lines: List[str],
) -> Tuple[List[_FeatureEntry], Dict[str, str]]:
    """
    Robust single-pass parser for the FEATURES section.

    Returns (list of _FeatureEntry, source_qualifiers_dict).
    """
    in_features = False
    features: List[_FeatureEntry] = []
    cur_feat: Optional[_FeatureEntry] = None
    source_qual: Dict[str, str] = {}

    # Accumulate multi-line location
    loc_buffer: List[str] = []
    in_loc = False

    for line in lines:
        if line.startswith('FEATURES'):
            in_features = True
            continue
        if not in_features:
            continue

        # Feature key line: 5 spaces + key + 2+ spaces + location
        feat_match = re.match(r'^     (\S+)\s{2,}(.+)$', line)
        if feat_match:
            # Finish previous feature
            if cur_feat:
                if loc_buffer:
                    cur_feat.location += ''.join(loc_buffer)
                    loc_buffer = []
                cur_feat.finalise_qualifier()
                features.append(cur_feat)
            cur_feat = _FeatureEntry(feat_match.group(1),
                                     feat_match.group(2).strip())
            in_loc = not _loc_complete(cur_feat.location)
            continue

        if cur_feat is None:
            continue

        # Location continuation or qualifier
        stripped = line.strip()

        if in_loc:
            # Still accumulating location
            cur_feat.location += stripped
            if _loc_complete(cur_feat.location):
                in_loc = False
            continue

        # Qualifier line: /key or /key="value"
        qual_match = re.match(r'^/(\w+)(=(.*))?$', stripped)
        if qual_match:
            cur_feat.finalise_qualifier()
            qkey = qual_match.group(1)
            qval = qual_match.group(3) or ''
            qval = qval.lstrip('"')
            if qval.endswith('"'):
                qval = qval[:-1]
                cur_feat._cur_qual = qkey
                cur_feat._cur_val  = [qval]
                cur_feat.finalise_qualifier()
            else:
                cur_feat._cur_qual = qkey
                cur_feat._cur_val  = [qval]
        elif stripped and cur_feat._cur_qual:
            # Continuation of a multi-line qualifier value
            chunk = stripped
            if chunk.endswith('"'):
                chunk = chunk[:-1]
                cur_feat._cur_val.append(' ' + chunk)
                cur_feat.finalise_qualifier()
            else:
                cur_feat._cur_val.append(' ' + chunk)

    # Append last feature
    if cur_feat:
        cur_feat.finalise_qualifier()
        features.append(cur_feat)

    # Extract source qualifiers
    for feat in features:
        if feat.key == 'source':
            source_qual = feat.qualifiers
            break

    return features, source_qual


def _loc_complete(loc: str) -> bool:
    """Return True if brackets are balanced (location string is complete)."""
    return loc.count('(') == loc.count(')')


def _process_feature_line(line, features, cur_feat):
    """Stub – replaced by _parse_features_singlepass."""
    pass


def _build_gene_models(
    features: List[_FeatureEntry],
    seq_id: str,
) -> List[GeneModel]:
    """
    Build GeneModel objects from a list of parsed _FeatureEntry objects.

    Strategy:
    1. Collect all gene, mRNA, and CDS features.
    2. Group mRNA and CDS by locus_tag or gene name.
    3. Assemble GeneModel for each gene.
    """
    gene_feats: Dict[str, _FeatureEntry] = {}
    mrna_feats: Dict[str, List[_FeatureEntry]] = {}
    cds_feats:  Dict[str, List[_FeatureEntry]] = {}

    def _gene_key(feat: _FeatureEntry) -> str:
        return (feat.qualifiers.get('locus_tag')
                or feat.qualifiers.get('gene')
                or feat.qualifiers.get('gene_synonym')
                or '')

    for feat in features:
        key = _gene_key(feat)
        if feat.key == 'gene':
            if key:
                gene_feats[key] = feat
        elif feat.key in ('mRNA', 'transcript'):
            mrna_feats.setdefault(key, []).append(feat)
        elif feat.key == 'CDS':
            cds_feats.setdefault(key, []).append(feat)

    models: List[GeneModel] = []

    # Process CDS features (the anchor for gene models)
    for gene_key, cds_list in cds_feats.items():
        for cds_idx, cds_feat in enumerate(cds_list):
            t_id = (cds_feat.qualifiers.get('transcript_id')
                    or f'{gene_key}.t{cds_idx + 1}')

            ivs, strand, p5, p3 = _parse_location(cds_feat.location)

            # Determine gene span
            if gene_key in gene_feats:
                g_ivs, _, _, _ = _parse_location(gene_feats[gene_key].location)
                gs = min(s for s, _ in g_ivs)
                ge = max(e for _, e in g_ivs)
            else:
                gs = min(s for s, _ in ivs) if ivs else 0
                ge = max(e for _, e in ivs) if ivs else 0

            cds_intervals = [CDSInterval(s, e, 0) for s, e in ivs]

            # Gather mRNA exons for this gene key
            exon_intervals: List[ExonInterval] = []
            if gene_key in mrna_feats:
                mrna = mrna_feats[gene_key][cds_idx] \
                    if cds_idx < len(mrna_feats[gene_key]) \
                    else mrna_feats[gene_key][0]
                mrna_ivs, _, _, _ = _parse_location(mrna.location)
                exon_intervals = [ExonInterval(s, e) for s, e in mrna_ivs]

            # Qualifiers
            qualifiers: Dict[str, str] = {}
            for qk in ('gene', 'locus_tag', 'product', 'protein_id',
                        'note', 'pseudo'):
                val = cds_feat.qualifiers.get(qk)
                if val:
                    qualifiers[qk] = val

            # codon_start
            cs = cds_feat.qualifiers.get('codon_start', '1')
            try:
                codon_start = int(cs)
            except ValueError:
                codon_start = 1
            if codon_start != 1:
                qualifiers['codon_start'] = str(codon_start)

            # translation → protein_seq
            translation = cds_feat.qualifiers.get('translation', '')
            # translation may be stored with internal spaces added during parsing
            translation = translation.replace(' ', '')

            gene_id = gene_key or cds_feat.qualifiers.get('gene', f'gene_{len(models)+1}')

            model = GeneModel(
                gene_id        = gene_id,
                transcript_id  = t_id,
                seq_id         = seq_id,
                strand         = strand,
                gene_start     = gs,
                gene_stop      = ge,
                cds_intervals  = cds_intervals,
                exon_intervals = exon_intervals,
                is_partial_5   = p5,
                is_partial_3   = p3,
                protein_seq    = translation or None,
                qualifiers     = qualifiers,
            )
            models.append(model)

    return models
