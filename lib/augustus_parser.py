"""
augustus_parser.py
------------------
Parse AUGUSTUS web-server HTML output or plain-text AUGUSTUS GFF output.

Supported input formats
-----------------------
1. AUGUSTUS HTML (output from https://bioinf.uni-greifswald.de/augustus/):
   The HTML contains a <pre class="result"> block with GFF predictions
   followed by inline coding and protein sequences as comment lines.

2. Plain AUGUSTUS GFF text output (the same text that appears inside the
   HTML <pre> block, saved as a plain .txt or .gff file).

AUGUSTUS GFF format summary
----------------------------
Columns (tab-separated):
  seqname   source  feature  start  end  score  strand  frame  attributes

Features used:
  gene        - gene span
  transcript  - transcript span (= gene when single-isoform)
  CDS         - coding exon(s)
  exon        - all exons including UTR (if predicted)
  start_codon - marks complete 5' end
  stop_codon  - marks complete 3' end

Attribute field (column 9):
  transcript_id "g1.t1"; gene_id "g1";
  -or- just: g1.t1  (for transcript/gene lines)

Inline sequence comments (appear AFTER all GFF lines for a transcript):
  # coding sequence = [ATGCGT...]
  # CONTINUATION...]
  # protein sequence = [MASQNKEE...]
  # CONTINUATION...]
"""

import re
from html.parser import HTMLParser
from typing import Dict, List, Optional, Tuple

from .models import GeneModel, CDSInterval, ExonInterval


# ─────────────────────────────────────────────────────────────────────────────
# HTML stripper – extracts raw text from AUGUSTUS HTML
# ─────────────────────────────────────────────────────────────────────────────

class _PreExtractor(HTMLParser):
    """Extract the content of the first <pre class="result"> element."""

    def __init__(self):
        super().__init__()
        self._in_pre = False
        self._depth  = 0
        self.text    = []

    def handle_starttag(self, tag, attrs):
        if tag == 'pre':
            attr_dict = dict(attrs)
            if attr_dict.get('class', '') == 'result':
                self._in_pre = True
                self._depth  = 0
            elif self._in_pre:
                self._depth += 1

    def handle_endtag(self, tag):
        if tag == 'pre' and self._in_pre:
            if self._depth == 0:
                self._in_pre = False
            else:
                self._depth -= 1

    def handle_data(self, data):
        if self._in_pre:
            self.text.append(data)

    def handle_entityref(self, name):
        if self._in_pre:
            entities = {'lt': '<', 'gt': '>', 'amp': '&', 'quot': '"'}
            self.text.append(entities.get(name, ''))

    def handle_charref(self, name):
        if self._in_pre:
            try:
                if name.startswith('x'):
                    self.text.append(chr(int(name[1:], 16)))
                else:
                    self.text.append(chr(int(name)))
            except ValueError:
                pass


def _extract_pre_text(html_path: str) -> str:
    """
    Return the plain text inside <pre class="result"> from an AUGUSTUS
    HTML file.  Falls back to returning the entire file content if no
    such element is found (handles plain-text GFF files).
    """
    with open(html_path, 'r', encoding='utf-8', errors='replace') as fh:
        content = fh.read()

    extractor = _PreExtractor()
    extractor.feed(content)
    text = ''.join(extractor.text)

    if not text.strip():
        # No <pre class="result"> found – treat entire file as plain text GFF
        return content
    return text


# ─────────────────────────────────────────────────────────────────────────────
# GFF attribute parser
# ─────────────────────────────────────────────────────────────────────────────

_ATTR_KEY_VAL = re.compile(r'(\w+)\s+"([^"]+)"')
_ATTR_PLAIN   = re.compile(r'^(\S+)$')


def _parse_attributes(attr_str: str) -> Dict[str, str]:
    """
    Parse AUGUSTUS attribute column.
    Handles both:
      transcript_id "g1.t1"; gene_id "g1";
    and plain (for gene/transcript lines):
      g1
      g1.t1
    """
    attr = {}
    matches = _ATTR_KEY_VAL.findall(attr_str)
    if matches:
        for key, val in matches:
            attr[key] = val
    else:
        plain = attr_str.strip()
        if '.' in plain:
            attr['transcript_id'] = plain
            attr['gene_id'] = plain.rsplit('.', 1)[0]
        elif plain:
            attr['gene_id'] = plain
    return attr


# ─────────────────────────────────────────────────────────────────────────────
# Inline sequence parser
# ─────────────────────────────────────────────────────────────────────────────

def _extract_inline_sequence(lines: List[str], start_idx: int,
                              keyword: str) -> Tuple[str, int]:
    """
    Extract a multi-line inline sequence from AUGUSTUS comment lines.

    AUGUSTUS format:
        # coding sequence = [ATGCGT...
        # MOREBASES...]
    or single line:
        # coding sequence = [ATGCGT...]

    Parameters
    ----------
    lines     : list of strings (all lines in the GFF block)
    start_idx : index of the line containing the keyword
    keyword   : 'coding sequence' or 'protein sequence'

    Returns
    -------
    (sequence_string, next_line_index)
    """
    seq_parts = []
    line = lines[start_idx].strip()

    # Find the opening bracket
    bracket_pos = line.find('[')
    if bracket_pos == -1:
        return '', start_idx + 1

    content = line[bracket_pos + 1:]
    # Check if closed on same line
    close_pos = content.find(']')
    if close_pos != -1:
        seq_parts.append(content[:close_pos])
        return ''.join(seq_parts), start_idx + 1

    seq_parts.append(content)
    idx = start_idx + 1

    while idx < len(lines):
        ln = lines[idx].strip()
        if ln.startswith('#'):
            text = ln[1:].strip()
            close_pos = text.find(']')
            if close_pos != -1:
                seq_parts.append(text[:close_pos])
                return ''.join(seq_parts), idx + 1
            else:
                seq_parts.append(text)
            idx += 1
        else:
            break

    return ''.join(seq_parts), idx


# ─────────────────────────────────────────────────────────────────────────────
# Main parser
# ─────────────────────────────────────────────────────────────────────────────

def parse_augustus_output(filepath: str) -> List[GeneModel]:
    """
    Parse an AUGUSTUS HTML or plain-text GFF file and return a list of
    GeneModel objects.

    Processing steps
    ----------------
    1. Extract the GFF text (strip HTML if needed).
    2. Parse GFF lines to build a dictionary of transcript models keyed by
       transcript_id.
    3. Parse inline '# coding sequence' and '# protein sequence' comment
       blocks and attach them to the most recently seen transcript.
    4. Determine partiality: a gene is 5'-partial if no start_codon line
       is present; 3'-partial if no stop_codon line is present.
    5. Return one GeneModel per transcript.

    Notes
    -----
    When AUGUSTUS is run with ``--alternatives-from-sampling=true``, there
    can be multiple transcripts per gene (e.g. g1.t1 and g1.t2).  All
    isoforms are returned; use the ``--primary-only`` flag in the CLI tools
    to keep only the first isoform (t1) per gene.

    Coordinate convention
    ---------------------
    All coordinates in the returned GeneModel objects are 1-based and
    stored with start <= stop, regardless of strand.
    """
    raw_text = _extract_pre_text(filepath)
    lines = raw_text.splitlines()

    # ── Data structures ──────────────────────────────────────────────────
    # gene_spans[gene_id]       = (seq_id, start, stop, strand, score)
    # transcript_data[t_id]     = dict of accumulated data
    # transcript_order          = insertion-ordered list of transcript_ids

    gene_spans: Dict[str, tuple]  = {}
    transcript_data: Dict[str, dict] = {}
    transcript_order: List[str]   = []
    current_transcript_id: Optional[str] = None

    i = 0
    while i < len(lines):
        line = lines[i]

        # ── Skip comment lines (but parse inline sequences) ───────────────
        if line.startswith('#'):
            stripped = line.strip()

            if '# coding sequence' in stripped and '[' in stripped:
                seq, i = _extract_inline_sequence(lines, i, 'coding sequence')
                if current_transcript_id and current_transcript_id in transcript_data:
                    transcript_data[current_transcript_id]['coding_seq'] = seq
                continue

            if '# protein sequence' in stripped and '[' in stripped:
                seq, i = _extract_inline_sequence(lines, i, 'protein sequence')
                if current_transcript_id and current_transcript_id in transcript_data:
                    # Strip trailing '*' if present
                    transcript_data[current_transcript_id]['protein_seq'] = \
                        seq.rstrip('*')
                continue

            i += 1
            continue

        # ── Parse GFF data line ───────────────────────────────────────────
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 8:
            i += 1
            continue

        seq_id   = cols[0]
        feature  = cols[2]
        try:
            start    = int(cols[3])
            stop     = int(cols[4])
        except ValueError:
            i += 1
            continue
        score_str = cols[5]
        strand   = cols[6]
        frame_str = cols[7]
        attr_str = cols[8] if len(cols) > 8 else ''

        try:
            score = float(score_str) if score_str not in ('.', '') else 0.0
        except ValueError:
            score = 0.0

        try:
            frame = int(frame_str) if frame_str not in ('.', '') else 0
        except ValueError:
            frame = 0

        attrs = _parse_attributes(attr_str)
        gene_id       = attrs.get('gene_id', '')
        transcript_id = attrs.get('transcript_id', gene_id)

        # ── Handle each feature type ──────────────────────────────────────

        if feature == 'gene':
            gene_id_here = gene_id or attr_str.strip()
            gene_spans[gene_id_here] = (seq_id, start, stop, strand, score)
            i += 1
            continue

        if feature in ('transcript', 'mRNA'):
            t_id = transcript_id or gene_id
            if t_id not in transcript_data:
                transcript_data[t_id] = {
                    'seq_id':          seq_id,
                    'strand':          strand,
                    'gene_id':         gene_id,
                    'transcript_id':   t_id,
                    'gene_start':      start,
                    'gene_stop':       stop,
                    'score':           score,
                    'cds_intervals':   [],
                    'exon_intervals':  [],
                    'has_start_codon': False,
                    'has_stop_codon':  False,
                    'coding_seq':      None,
                    'protein_seq':     None,
                }
                transcript_order.append(t_id)
            current_transcript_id = t_id
            i += 1
            continue

        # For CDS/exon/start_codon/stop_codon – look up transcript_id
        if feature in ('CDS', 'exon', 'start_codon', 'stop_codon'):
            t_id = transcript_id
            # Create transcript entry if we haven't seen the transcript line
            # (this can happen in some GFF variants)
            if t_id and t_id not in transcript_data:
                transcript_data[t_id] = {
                    'seq_id':          seq_id,
                    'strand':          strand,
                    'gene_id':         gene_id,
                    'transcript_id':   t_id,
                    'gene_start':      start,
                    'gene_stop':       stop,
                    'score':           score,
                    'cds_intervals':   [],
                    'exon_intervals':  [],
                    'has_start_codon': False,
                    'has_stop_codon':  False,
                    'coding_seq':      None,
                    'protein_seq':     None,
                }
                transcript_order.append(t_id)
                current_transcript_id = t_id
            elif t_id:
                current_transcript_id = t_id

            if not t_id or t_id not in transcript_data:
                i += 1
                continue

            td = transcript_data[t_id]

            if feature == 'CDS':
                td['cds_intervals'].append(CDSInterval(start, stop, frame))

            elif feature == 'exon':
                td['exon_intervals'].append(ExonInterval(start, stop))

            elif feature == 'start_codon':
                td['has_start_codon'] = True

            elif feature == 'stop_codon':
                td['has_stop_codon'] = True

        i += 1

    # ── Build GeneModel objects ───────────────────────────────────────────
    gene_models: List[GeneModel] = []

    for t_id in transcript_order:
        td = transcript_data[t_id]

        if not td['cds_intervals']:
            # Skip transcripts with no CDS (non-coding or parse error)
            continue

        gene_id = td['gene_id']
        # Use gene span if available, else use transcript span
        if gene_id in gene_spans:
            _, gs, ge, _, _ = gene_spans[gene_id]
        else:
            gs = min(iv.start for iv in td['cds_intervals'])
            ge = max(iv.stop  for iv in td['cds_intervals'])

        is_partial_5 = not td['has_start_codon']
        is_partial_3 = not td['has_stop_codon']

        # codon_start qualifier for partial 5' CDS
        qualifiers: Dict[str, str] = {}
        if is_partial_5:
            # Determine frame of first CDS interval in transcript order
            sorted_cds = sorted(td['cds_intervals'], key=lambda iv: iv.start)
            if td['strand'] == '-':
                first_iv = sorted_cds[-1]
            else:
                first_iv = sorted_cds[0]
            phase = first_iv.phase
            # GFF phase → codon_start: phase 0→1, phase 1→2, phase 2→3
            codon_start = phase + 1 if phase in (0, 1, 2) else 1
            if codon_start > 1:
                qualifiers['codon_start'] = str(codon_start)

        model = GeneModel(
            gene_id        = gene_id,
            transcript_id  = t_id,
            seq_id         = td['seq_id'],
            strand         = td['strand'],
            gene_start     = gs,
            gene_stop      = ge,
            cds_intervals  = td['cds_intervals'],
            exon_intervals = td['exon_intervals'],
            is_partial_5   = is_partial_5,
            is_partial_3   = is_partial_3,
            coding_seq     = td['coding_seq'],
            protein_seq    = td['protein_seq'],
            score          = td['score'],
            qualifiers     = qualifiers,
        )
        gene_models.append(model)

    return gene_models


def filter_primary_isoforms(gene_models: List[GeneModel]) -> List[GeneModel]:
    """
    Keep only the first isoform (t1) for each gene_id.
    This is used when ``--primary-only`` is specified.
    """
    seen_genes: set = set()
    primary: List[GeneModel] = []
    for gm in gene_models:
        if gm.gene_id not in seen_genes:
            seen_genes.add(gm.gene_id)
            primary.append(gm)
    return primary
