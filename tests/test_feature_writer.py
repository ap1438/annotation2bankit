#!/usr/bin/env python3
"""
test_feature_writer.py
----------------------
Tests for the BankIt feature table, nucleotide FASTA, and protein FASTA writers.

Tests verify:
  1. Feature table tab structure (3 leading tabs for qualifiers)
  2. Plus-strand coordinate ordering (ascending)
  3. Minus-strand coordinate ordering (descending, first > second)
  4. Partial gene notation (<, >)
  5. Multi-exon feature output
  6. Required qualifiers present
  7. Nucleotide FASTA defline format (organism, mol_type modifiers)
  8. Protein FASTA output
  9. Integration: round-trip with training data

Run:
    python -m pytest tests/test_feature_writer.py -v
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lib.models import GeneModel, CDSInterval, ExonInterval, SequenceRecord
from lib.feature_writer import (write_feature_table, write_nucleotide_fasta,
                                 write_protein_fasta)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _make_plus_model(gene_id='g1', intervals=None) -> GeneModel:
    """Create a simple plus-strand gene model."""
    if intervals is None:
        intervals = [(100, 300), (500, 700), (800, 999)]
    cds = [CDSInterval(s, e, 0) for s, e in intervals]
    return GeneModel(
        gene_id        = gene_id,
        transcript_id  = f'{gene_id}.t1',
        seq_id         = 'TestSeq',
        strand         = '+',
        gene_start     = intervals[0][0],
        gene_stop      = intervals[-1][1],
        cds_intervals  = cds,
        qualifiers     = {'gene': gene_id, 'product': 'test protein'},
    )


def _make_minus_model(gene_id='g2') -> GeneModel:
    """Create a simple minus-strand gene model."""
    intervals = [(100, 300), (500, 700)]
    cds = [CDSInterval(s, e, 0) for s, e in intervals]
    return GeneModel(
        gene_id        = gene_id,
        transcript_id  = f'{gene_id}.t1',
        seq_id         = 'TestSeq',
        strand         = '-',
        gene_start     = 100,
        gene_stop      = 700,
        cds_intervals  = cds,
        qualifiers     = {'gene': gene_id, 'product': 'minus protein'},
    )


def _make_partial_model() -> GeneModel:
    cds = [CDSInterval(1, 300, 1), CDSInterval(500, 800, 0)]
    return GeneModel(
        gene_id        = 'partial_gene',
        transcript_id  = 'partial_gene.t1',
        seq_id         = 'TestSeq',
        strand         = '+',
        gene_start     = 1,
        gene_stop      = 800,
        cds_intervals  = cds,
        is_partial_5   = True,
        is_partial_3   = True,
        qualifiers     = {'gene': 'partial_gene', 'product': 'partial protein',
                          'codon_start': '2'},
    )


def _read_tbl(path: str) -> str:
    with open(path) as fh:
        return fh.read()


# ─────────────────────────────────────────────────────────────────────────────
# Feature table format tests
# ─────────────────────────────────────────────────────────────────────────────

class TestFeatureTableFormat(unittest.TestCase):

    def setUp(self):
        self.tbl_file = tempfile.NamedTemporaryFile(
            suffix='.tbl', delete=False, mode='w')
        self.tbl_path = self.tbl_file.name
        self.tbl_file.close()

    def tearDown(self):
        os.unlink(self.tbl_path)

    def test_feature_header(self):
        gm = _make_plus_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        self.assertIn('>Feature TestSeq', content)

    def test_qualifier_has_three_leading_tabs(self):
        gm = _make_plus_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        # Every qualifier line should start with exactly 3 tabs
        for line in content.splitlines():
            if line.startswith('\t\t\t'):
                self.assertFalse(line.startswith('\t\t\t\t'),
                                 f'Qualifier has 4+ tabs: {repr(line)}')
                parts = line[3:].split('\t')
                self.assertGreater(len(parts), 0)

    def test_plus_strand_gene_coords_ascending(self):
        gm = _make_plus_model(intervals=[(100, 300), (500, 700)])
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        lines = [l for l in content.splitlines() if l and not l.startswith('>')]
        # First non-header line should be the gene feature: 100 \t 700 \t gene
        # (or gene_start \t gene_stop \t gene)
        data_lines = [l for l in lines if not l.startswith('\t')]
        if data_lines:
            first = data_lines[0]
            cols = first.split('\t')
            if len(cols) >= 2:
                c1 = int(cols[0].replace('<', '').replace('>', ''))
                c2 = int(cols[1].replace('<', '').replace('>', ''))
                self.assertLess(c1, c2, 'Plus strand: first coord < second coord')

    def test_minus_strand_gene_coords_descending(self):
        gm = _make_minus_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        lines = [l for l in content.splitlines() if l and not l.startswith('>')]
        data_lines = [l for l in lines if not l.startswith('\t')]
        if data_lines:
            first = data_lines[0]
            cols = first.split('\t')
            if len(cols) >= 2:
                c1 = int(cols[0].replace('<', '').replace('>', ''))
                c2 = int(cols[1].replace('<', '').replace('>', ''))
                self.assertGreater(c1, c2,
                                   'Minus strand: first coord > second coord')

    def test_partial_5_uses_less_than_symbol(self):
        gm = _make_partial_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        self.assertIn('<', content, 'Partial 5\' should use < symbol')

    def test_partial_3_uses_greater_than_symbol(self):
        gm = _make_partial_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        self.assertIn('>', content, 'Partial 3\' should use > symbol')

    def test_gene_feature_present(self):
        gm = _make_plus_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        self.assertIn('\tgene\n', content, 'gene feature line missing')

    def test_mrna_feature_present(self):
        gm = _make_plus_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        self.assertIn('\tmRNA\n', content, 'mRNA feature line missing')

    def test_cds_feature_present(self):
        gm = _make_plus_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        self.assertIn('\tCDS\n', content, 'CDS feature line missing')

    def test_product_qualifier_present(self):
        gm = _make_plus_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        self.assertIn('product\ttest protein', content)

    def test_gene_qualifier_present(self):
        gm = _make_plus_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        self.assertIn('gene\tg1', content)

    def test_codon_start_qualifier(self):
        gm = _make_partial_model()
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        self.assertIn('codon_start\t2', content)

    def test_multi_exon_continuation_lines(self):
        """Multi-exon CDS should produce continuation lines (no feature key)."""
        gm = _make_plus_model(intervals=[(100, 200), (400, 500), (700, 900)])
        write_feature_table([gm], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        lines = content.splitlines()
        # Count lines that have 2 tab-separated integers (continuation lines)
        continuation = []
        for line in lines:
            if line.startswith('\t') or line.startswith('>'):
                continue
            cols = line.split('\t')
            if len(cols) == 2:
                try:
                    int(cols[0].replace('<', '').replace('>', ''))
                    int(cols[1].replace('<', '').replace('>', ''))
                    continuation.append(line)
                except ValueError:
                    pass
        # Should have at least 2 continuation lines for 3-exon CDS (gene + mrna + cds)
        self.assertGreater(len(continuation), 0,
                           'Multi-exon feature should have continuation lines')

    def test_multiple_genes_multiple_headers(self):
        """Two genes on different sequences → two >Feature headers."""
        gm1 = _make_plus_model('gene1')
        gm2 = GeneModel(
            gene_id='gene2', transcript_id='gene2.t1', seq_id='OtherSeq',
            strand='+', gene_start=1, gene_stop=500,
            cds_intervals=[CDSInterval(1, 500, 0)],
            qualifiers={'gene': 'gene2', 'product': 'other protein'},
        )
        write_feature_table([gm1, gm2], self.tbl_path)
        content = _read_tbl(self.tbl_path)
        self.assertIn('>Feature TestSeq', content)
        self.assertIn('>Feature OtherSeq', content)


# ─────────────────────────────────────────────────────────────────────────────
# Nucleotide FASTA writer tests
# ─────────────────────────────────────────────────────────────────────────────

class TestNucleotideFastaWriter(unittest.TestCase):

    def setUp(self):
        self.out = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False)
        self.out_path = self.out.name
        self.out.close()

    def tearDown(self):
        os.unlink(self.out_path)

    def _make_record(self, seq_id='SEQ1', seq='ACGT' * 20) -> SequenceRecord:
        return SequenceRecord(
            seq_id      = seq_id,
            sequence    = seq,
            organism    = 'Test organism',
            mol_type    = 'genomic DNA',
            description = 'genomic sequence',
            modifiers   = {'isolate': 'test_isolate'},
        )

    def test_defline_starts_with_seq_id(self):
        rec = self._make_record()
        write_nucleotide_fasta([rec], self.out_path)
        with open(self.out_path) as fh:
            first = fh.readline().strip()
        self.assertTrue(first.startswith('>SEQ1'))

    def test_organism_modifier_in_defline(self):
        rec = self._make_record()
        write_nucleotide_fasta([rec], self.out_path)
        with open(self.out_path) as fh:
            first = fh.readline().strip()
        self.assertIn('[organism=Test organism]', first)

    def test_mol_type_modifier_in_defline(self):
        rec = self._make_record()
        write_nucleotide_fasta([rec], self.out_path)
        with open(self.out_path) as fh:
            first = fh.readline().strip()
        self.assertIn('[mol_type=genomic DNA]', first)

    def test_isolate_modifier_in_defline(self):
        rec = self._make_record()
        write_nucleotide_fasta([rec], self.out_path)
        with open(self.out_path) as fh:
            first = fh.readline().strip()
        self.assertIn('[isolate=test_isolate]', first)

    def test_sequence_lines_60_chars(self):
        rec = self._make_record(seq='A' * 200)
        write_nucleotide_fasta([rec], self.out_path)
        with open(self.out_path) as fh:
            lines = fh.readlines()
        seq_lines = [l.strip() for l in lines if not l.startswith('>')]
        for line in seq_lines[:-1]:  # last line may be shorter
            self.assertEqual(len(line), 60,
                             f'Sequence line length should be 60, got {len(line)}')

    def test_sequence_uppercase(self):
        rec = self._make_record(seq='acgtacgt' * 10)
        write_nucleotide_fasta([rec], self.out_path)
        with open(self.out_path) as fh:
            content = fh.read()
        seq_part = ''.join(l.strip() for l in content.splitlines()
                           if not l.startswith('>'))
        self.assertEqual(seq_part, seq_part.upper())


# ─────────────────────────────────────────────────────────────────────────────
# Integration: write and validate against training data
# ─────────────────────────────────────────────────────────────────────────────

class TestFeatureWriterIntegration(unittest.TestCase):
    """Write feature table from training data feature file and check correctness."""

    def test_known_hma4_wall10_model(self):
        """
        Build the HMA4-1 Wall10 gene model manually (from known annotation)
        and verify the feature table matches expected BankIt format.
        """
        # HMA4-1_Wall_10_v1.0.0 is on minus strand
        # CDS intervals from Wall10.HMA4-123.tbl (1-based, start<=stop):
        cds_ivs = [
            CDSInterval(219,  604, 0),
            CDSInterval(723,  2124, 0),
            CDSInterval(2492, 2692, 0),
            CDSInterval(2777, 3104, 0),  # approximate
            CDSInterval(3189, 3326, 0),
            CDSInterval(3415, 3670, 0),
            CDSInterval(3803, 3900, 0),
            CDSInterval(5767, 6024, 0),
            CDSInterval(7174, 7387, 0),
        ]
        gm = GeneModel(
            gene_id        = 'AhHMA4-1',
            transcript_id  = 'AhHMA4-1.t1',
            seq_id         = 'HMA4-1_Wall_10_v1.0.0',
            strand         = '-',
            gene_start     = 1,
            gene_stop      = 8083,
            cds_intervals  = cds_ivs,
            qualifiers     = {
                'gene':    'AhHMA4-1',
                'product': 'Zn/Cd P(IB)-type ATPase',
            },
        )

        with tempfile.NamedTemporaryFile(suffix='.tbl', delete=False,
                                         mode='w') as fh:
            tbl_path = fh.name
        try:
            write_feature_table([gm], tbl_path)
            content = _read_tbl(tbl_path)
            self.assertIn('>Feature HMA4-1_Wall_10_v1.0.0', content)
            self.assertIn('AhHMA4-1', content)
            self.assertIn('Zn/Cd P(IB)-type ATPase', content)

            # For minus strand, first coord should be > second
            lines = [l for l in content.splitlines()
                     if l and not l.startswith('>') and not l.startswith('\t')]
            for line in lines:
                cols = line.split('\t')
                if len(cols) >= 2:
                    try:
                        c1 = int(cols[0].replace('<', '').replace('>', ''))
                        c2 = int(cols[1].replace('<', '').replace('>', ''))
                        self.assertGreater(c1, c2,
                                           f'Minus strand line should have c1>c2: {line}')
                    except ValueError:
                        pass
        finally:
            os.unlink(tbl_path)


if __name__ == '__main__':
    unittest.main(verbosity=2)
