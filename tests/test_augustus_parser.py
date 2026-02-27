#!/usr/bin/env python3
"""
test_augustus_parser.py
-----------------------
Unit and integration tests for the AUGUSTUS HTML/GFF parser.

Tests verify:
  1. AUGUSTUS HTML parsing produces correct gene models
  2. CDS intervals are correctly extracted
  3. Strand is correctly detected
  4. Coding and protein sequences are extracted
  5. Partiality flags are set correctly
  6. Primary isoform filtering works

Run:
    python -m pytest tests/test_augustus_parser.py -v
    # or
    python tests/test_augustus_parser.py
"""

import os
import sys
import tempfile
import unittest

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lib.augustus_parser import parse_augustus_output, filter_primary_isoforms
from lib.models import GeneModel, CDSInterval


# ─────────────────────────────────────────────────────────────────────────────
# Minimal AUGUSTUS GFF for testing
# ─────────────────────────────────────────────────────────────────────────────

PLUS_STRAND_GFF = """\
# AUGUSTUS prediction
# ----- prediction on sequence number 1 (length = 10000, name = TestSeq) -----
TestSeq\tAUGUSTUS\tgene\t100\t900\t0.9\t+\t.\tg1
TestSeq\tAUGUSTUS\ttranscript\t100\t900\t.\t+\t.\tg1.t1
TestSeq\tAUGUSTUS\tstart_codon\t100\t102\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq\tAUGUSTUS\tCDS\t100\t300\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq\tAUGUSTUS\texon\t100\t300\t.\t+\t.\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq\tAUGUSTUS\tCDS\t400\t600\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq\tAUGUSTUS\texon\t400\t600\t.\t+\t.\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq\tAUGUSTUS\tCDS\t700\t900\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq\tAUGUSTUS\texon\t700\t900\t.\t+\t.\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq\tAUGUSTUS\tstop_codon\t898\t900\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
# coding sequence = [ATGCGT]
# protein sequence = [MR]
"""

MINUS_STRAND_GFF = """\
# AUGUSTUS prediction
TestSeq2\tAUGUSTUS\tgene\t1\t800\t0.8\t-\t.\tg1
TestSeq2\tAUGUSTUS\ttranscript\t1\t800\t.\t-\t.\tg1.t1
TestSeq2\tAUGUSTUS\tstart_codon\t598\t600\t.\t-\t0\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq2\tAUGUSTUS\tCDS\t400\t600\t.\t-\t0\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq2\tAUGUSTUS\tCDS\t100\t300\t.\t-\t0\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq2\tAUGUSTUS\tstop_codon\t100\t102\t.\t-\t0\ttranscript_id "g1.t1"; gene_id "g1";
# coding sequence = [ATGCGT]
# protein sequence = [MR]
"""

PARTIAL_GFF = """\
# AUGUSTUS prediction
TestSeq3\tAUGUSTUS\tgene\t1\t500\t0.7\t+\t.\tg1
TestSeq3\tAUGUSTUS\ttranscript\t1\t500\t.\t+\t.\tg1.t1
TestSeq3\tAUGUSTUS\tCDS\t1\t200\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";
TestSeq3\tAUGUSTUS\tCDS\t300\t500\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
# coding sequence = [CGTAAA]
# protein sequence = [RK]
"""

MULTI_GENE_GFF = """\
# AUGUSTUS prediction
Seq\tAUGUSTUS\tgene\t100\t500\t0.9\t+\t.\tg1
Seq\tAUGUSTUS\ttranscript\t100\t500\t.\t+\t.\tg1.t1
Seq\tAUGUSTUS\tstart_codon\t100\t102\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
Seq\tAUGUSTUS\tCDS\t100\t500\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
Seq\tAUGUSTUS\tstop_codon\t498\t500\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
# coding sequence = [ATGTAA]
# protein sequence = [M]
Seq\tAUGUSTUS\tgene\t1000\t2000\t0.85\t-\t.\tg2
Seq\tAUGUSTUS\ttranscript\t1000\t2000\t.\t-\t.\tg2.t1
Seq\tAUGUSTUS\tstart_codon\t1998\t2000\t.\t-\t0\ttranscript_id "g2.t1"; gene_id "g2";
Seq\tAUGUSTUS\tCDS\t1500\t2000\t.\t-\t0\ttranscript_id "g2.t1"; gene_id "g2";
Seq\tAUGUSTUS\tCDS\t1000\t1300\t.\t-\t0\ttranscript_id "g2.t1"; gene_id "g2";
Seq\tAUGUSTUS\tstop_codon\t1000\t1002\t.\t-\t0\ttranscript_id "g2.t1"; gene_id "g2";
# coding sequence = [ATGGGG]
# protein sequence = [MG]
"""

MULTI_ISOFORM_GFF = """\
# AUGUSTUS prediction
Seq\tAUGUSTUS\tgene\t100\t900\t0.9\t+\t.\tg1
Seq\tAUGUSTUS\ttranscript\t100\t900\t.\t+\t.\tg1.t1
Seq\tAUGUSTUS\tstart_codon\t100\t102\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
Seq\tAUGUSTUS\tCDS\t100\t500\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
Seq\tAUGUSTUS\tstop_codon\t498\t500\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
# protein sequence = [MAAAA]
Seq\tAUGUSTUS\ttranscript\t100\t900\t.\t+\t.\tg1.t2
Seq\tAUGUSTUS\tstart_codon\t100\t102\t.\t+\t0\ttranscript_id "g1.t2"; gene_id "g1";
Seq\tAUGUSTUS\tCDS\t100\t300\t.\t+\t0\ttranscript_id "g1.t2"; gene_id "g1";
Seq\tAUGUSTUS\tCDS\t600\t900\t.\t+\t0\ttranscript_id "g1.t2"; gene_id "g1";
Seq\tAUGUSTUS\tstop_codon\t898\t900\t.\t+\t0\ttranscript_id "g1.t2"; gene_id "g1";
# protein sequence = [MAAAB]
"""


def _write_temp_gff(content: str) -> str:
    """Write content to a temp file and return path."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gff',
                                     delete=False) as fh:
        fh.write(content)
        return fh.name


class TestAugustusParserBasic(unittest.TestCase):
    """Basic parsing tests for plus-strand, minus-strand, and partial genes."""

    def test_plus_strand_gene(self):
        path = _write_temp_gff(PLUS_STRAND_GFF)
        try:
            models = parse_augustus_output(path)
            self.assertEqual(len(models), 1, 'Expected 1 gene model')
            gm = models[0]
            self.assertEqual(gm.gene_id, 'g1')
            self.assertEqual(gm.strand, '+')
            self.assertEqual(gm.seq_id, 'TestSeq')
            self.assertEqual(len(gm.cds_intervals), 3)
            self.assertFalse(gm.is_partial_5)
            self.assertFalse(gm.is_partial_3)
        finally:
            os.unlink(path)

    def test_minus_strand_gene(self):
        path = _write_temp_gff(MINUS_STRAND_GFF)
        try:
            models = parse_augustus_output(path)
            self.assertEqual(len(models), 1)
            gm = models[0]
            self.assertEqual(gm.strand, '-')
            self.assertEqual(len(gm.cds_intervals), 2)
            self.assertFalse(gm.is_partial_5)
            self.assertFalse(gm.is_partial_3)
        finally:
            os.unlink(path)

    def test_partial_gene(self):
        path = _write_temp_gff(PARTIAL_GFF)
        try:
            models = parse_augustus_output(path)
            self.assertEqual(len(models), 1)
            gm = models[0]
            self.assertTrue(gm.is_partial_5, '5-partial flag should be set')
            self.assertTrue(gm.is_partial_3, '3-partial flag should be set')
        finally:
            os.unlink(path)

    def test_inline_sequences_extracted(self):
        path = _write_temp_gff(PLUS_STRAND_GFF)
        try:
            models = parse_augustus_output(path)
            gm = models[0]
            self.assertIsNotNone(gm.coding_seq,  'coding_seq should be set')
            self.assertIsNotNone(gm.protein_seq, 'protein_seq should be set')
            self.assertIn('ATG', gm.coding_seq.upper())
            self.assertIn('M',   gm.protein_seq.upper())
        finally:
            os.unlink(path)

    def test_exon_intervals_extracted(self):
        path = _write_temp_gff(PLUS_STRAND_GFF)
        try:
            models = parse_augustus_output(path)
            gm = models[0]
            self.assertEqual(len(gm.exon_intervals), 3)
        finally:
            os.unlink(path)


class TestAugustusParserMultiGene(unittest.TestCase):
    """Tests for multiple genes and isoforms."""

    def test_multi_gene_count(self):
        path = _write_temp_gff(MULTI_GENE_GFF)
        try:
            models = parse_augustus_output(path)
            self.assertEqual(len(models), 2, 'Expected 2 gene models')
        finally:
            os.unlink(path)

    def test_multi_gene_strands(self):
        path = _write_temp_gff(MULTI_GENE_GFF)
        try:
            models = parse_augustus_output(path)
            strands = {gm.gene_id: gm.strand for gm in models}
            self.assertEqual(strands.get('g1'), '+')
            self.assertEqual(strands.get('g2'), '-')
        finally:
            os.unlink(path)

    def test_multi_isoform(self):
        path = _write_temp_gff(MULTI_ISOFORM_GFF)
        try:
            models = parse_augustus_output(path)
            self.assertEqual(len(models), 2, 'Both isoforms should be parsed')
            t_ids = [gm.transcript_id for gm in models]
            self.assertIn('g1.t1', t_ids)
            self.assertIn('g1.t2', t_ids)
        finally:
            os.unlink(path)

    def test_primary_isoform_filter(self):
        path = _write_temp_gff(MULTI_ISOFORM_GFF)
        try:
            models  = parse_augustus_output(path)
            primary = filter_primary_isoforms(models)
            self.assertEqual(len(primary), 1)
            self.assertEqual(primary[0].transcript_id, 'g1.t1')
        finally:
            os.unlink(path)


class TestAugustusCoordinateSorting(unittest.TestCase):
    """Test that CDS intervals are sorted correctly for BankIt output."""

    def test_plus_strand_cds_ascending(self):
        path = _write_temp_gff(PLUS_STRAND_GFF)
        try:
            models = parse_augustus_output(path)
            gm = models[0]
            sorted_ivs = gm.cds_sorted_for_table()
            starts = [iv.start for iv in sorted_ivs]
            self.assertEqual(starts, sorted(starts),
                             'Plus strand CDS should be ascending')
        finally:
            os.unlink(path)

    def test_minus_strand_cds_descending(self):
        path = _write_temp_gff(MINUS_STRAND_GFF)
        try:
            models = parse_augustus_output(path)
            gm = models[0]
            sorted_ivs = gm.cds_sorted_for_table()
            stops = [iv.stop for iv in sorted_ivs]
            self.assertEqual(stops, sorted(stops, reverse=True),
                             'Minus strand CDS should be in descending order')
        finally:
            os.unlink(path)


class TestLiveDataFile(unittest.TestCase):
    """
    Integration test against the actual AUGUSTUS HTML file in the
    Training_data directory.  Skipped if file not present.
    """

    TRAINING_HTML = (
        '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
        'Augustus_ to Feature_file_to_genbank_file/Training_data/'
        'Augustus_ result display.html'
    )
    WALL10_HTML = (
        '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
        'Augustus_ to Feature_file_to_genbank_file/Augustus_ Wall10123.html'
    )

    @unittest.skipUnless(
        os.path.isfile(
            '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
            'Augustus_ to Feature_file_to_genbank_file/Training_data/'
            'Augustus_ result display.html'),
        'Training data not found'
    )
    def test_bac7c17_html_parses(self):
        models = parse_augustus_output(self.TRAINING_HTML)
        self.assertGreater(len(models), 0, 'Should parse at least one gene')
        print(f'\n  [INFO] BAC_7C17 HTML: {len(models)} gene model(s) parsed.')

    @unittest.skipUnless(
        os.path.isfile(
            '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
            'Augustus_ to Feature_file_to_genbank_file/Augustus_ Wall10123.html'),
        'Wall10 HTML not found'
    )
    def test_wall10_html_parses(self):
        models = parse_augustus_output(self.WALL10_HTML)
        self.assertGreater(len(models), 0)
        print(f'\n  [INFO] Wall10 HTML: {len(models)} gene model(s) parsed.')
        seq_ids = set(gm.seq_id for gm in models)
        print(f'  [INFO] Sequence IDs: {seq_ids}')


if __name__ == '__main__':
    unittest.main(verbosity=2)
