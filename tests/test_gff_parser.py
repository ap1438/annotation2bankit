#!/usr/bin/env python3
"""
test_gff_parser.py
------------------
Tests for GFF3, GTF, and GFF2 parsers.

Tests:
  1. GFF3 format detection
  2. GTF format detection
  3. GFF3 gene→mRNA→CDS hierarchy
  4. GFF3 partial features
  5. GFF3 minus strand
  6. GTF parsing
  7. GFF2 fallback parsing
  8. Embedded FASTA in GFF3
  9. Integration: BAC_7C17 GFF3

Run:
    python -m pytest tests/test_gff_parser.py -v
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lib.gff_parser import parse_gff, detect_format


# ─────────────────────────────────────────────────────────────────────────────
# Test GFF content
# ─────────────────────────────────────────────────────────────────────────────

GFF3_BASIC = """\
##gff-version 3
##sequence-region Chr1 1 10000
Chr1\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene1;Name=TestGene;locus_tag=TG001
Chr1\tRefSeq\tmRNA\t100\t900\t.\t+\t.\tID=mrna1;Parent=gene1;product=test protein
Chr1\tRefSeq\texon\t100\t300\t.\t+\t.\tID=exon1;Parent=mrna1
Chr1\tRefSeq\texon\t500\t700\t.\t+\t.\tID=exon2;Parent=mrna1
Chr1\tRefSeq\texon\t800\t900\t.\t+\t.\tID=exon3;Parent=mrna1
Chr1\tRefSeq\tCDS\t100\t300\t.\t+\t0\tID=cds1;Parent=mrna1
Chr1\tRefSeq\tCDS\t500\t700\t.\t+\t0\tID=cds2;Parent=mrna1
Chr1\tRefSeq\tCDS\t800\t900\t.\t+\t0\tID=cds3;Parent=mrna1
"""

GFF3_MINUS = """\
##gff-version 3
Chr2\tRefSeq\tgene\t100\t900\t.\t-\t.\tID=gene2;Name=MinusGene
Chr2\tRefSeq\tmRNA\t100\t900\t.\t-\t.\tID=mrna2;Parent=gene2;product=minus strand protein
Chr2\tRefSeq\tCDS\t500\t900\t.\t-\t0\tID=cds-a;Parent=mrna2
Chr2\tRefSeq\tCDS\t100\t300\t.\t-\t0\tID=cds-b;Parent=mrna2
"""

GFF3_PARTIAL = """\
##gff-version 3
Chr3\tRefSeq\tgene\t1\t500\t.\t+\t.\tID=gene3;Name=PartialGene;partial=true
Chr3\tRefSeq\tmRNA\t1\t500\t.\t+\t.\tID=mrna3;Parent=gene3;partial=true;start_range=.,1;product=partial protein
Chr3\tRefSeq\tCDS\t1\t500\t.\t+\t1\tID=cds3;Parent=mrna3
"""

GFF3_WITH_FASTA = """\
##gff-version 3
TestSeq\tRefSeq\tgene\t1\t300\t.\t+\t.\tID=g1;Name=FaGene
TestSeq\tRefSeq\tmRNA\t1\t300\t.\t+\t.\tID=m1;Parent=g1;product=fasta protein
TestSeq\tRefSeq\tCDS\t1\t300\t.\t+\t0\tID=c1;Parent=m1
##FASTA
>TestSeq
ATGCGTAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAAT
TTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAA
TTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAA
ATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGG
ATGA
"""

GTF_BASIC = """\
# GTF file
Chr1\tRefSeq\tgene\t100\t900\t.\t+\t.\tgene_id "GTFGene1";
Chr1\tRefSeq\ttranscript\t100\t900\t.\t+\t.\tgene_id "GTFGene1"; transcript_id "GTFGene1.1"; product "gtf protein";
Chr1\tRefSeq\texon\t100\t300\t.\t+\t.\tgene_id "GTFGene1"; transcript_id "GTFGene1.1";
Chr1\tRefSeq\texon\t500\t900\t.\t+\t.\tgene_id "GTFGene1"; transcript_id "GTFGene1.1";
Chr1\tRefSeq\tCDS\t100\t300\t.\t+\t0\tgene_id "GTFGene1"; transcript_id "GTFGene1.1"; product "gtf protein";
Chr1\tRefSeq\tCDS\t500\t900\t.\t+\t0\tgene_id "GTFGene1"; transcript_id "GTFGene1.1";
"""

GFF2_BASIC = """\
# GFF2 file
Chr1\tRefSeq\tCDS\t100\t300\t.\t+\t0\tgene_id "GFF2Gene1"
Chr1\tRefSeq\tCDS\t500\t700\t.\t+\t0\tgene_id "GFF2Gene1"
"""


def _write_temp(content: str, suffix: str = '.gff3') -> str:
    with tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False) as fh:
        fh.write(content)
        return fh.name


# ─────────────────────────────────────────────────────────────────────────────
# Format detection tests
# ─────────────────────────────────────────────────────────────────────────────

class TestFormatDetection(unittest.TestCase):

    def test_detect_gff3(self):
        path = _write_temp(GFF3_BASIC, '.gff3')
        try:
            self.assertEqual(detect_format(path), 'gff3')
        finally:
            os.unlink(path)

    def test_detect_gtf(self):
        path = _write_temp(GTF_BASIC, '.gtf')
        try:
            fmt = detect_format(path)
            self.assertEqual(fmt, 'gtf')
        finally:
            os.unlink(path)


# ─────────────────────────────────────────────────────────────────────────────
# GFF3 parsing tests
# ─────────────────────────────────────────────────────────────────────────────

class TestGFF3Parser(unittest.TestCase):

    def test_basic_parse(self):
        path = _write_temp(GFF3_BASIC)
        try:
            seqs, models = parse_gff(path, fmt='gff3')
            self.assertEqual(len(models), 1)
        finally:
            os.unlink(path)

    def test_gene_name(self):
        path = _write_temp(GFF3_BASIC)
        try:
            _, models = parse_gff(path, fmt='gff3')
            gm = models[0]
            self.assertEqual(gm.gene_id, 'TestGene')
        finally:
            os.unlink(path)

    def test_locus_tag(self):
        path = _write_temp(GFF3_BASIC)
        try:
            _, models = parse_gff(path, fmt='gff3')
            gm = models[0]
            self.assertEqual(gm.qualifiers.get('locus_tag'), 'TG001')
        finally:
            os.unlink(path)

    def test_product(self):
        path = _write_temp(GFF3_BASIC)
        try:
            _, models = parse_gff(path, fmt='gff3')
            gm = models[0]
            self.assertEqual(gm.qualifiers.get('product'), 'test protein')
        finally:
            os.unlink(path)

    def test_cds_interval_count(self):
        path = _write_temp(GFF3_BASIC)
        try:
            _, models = parse_gff(path, fmt='gff3')
            gm = models[0]
            self.assertEqual(len(gm.cds_intervals), 3)
        finally:
            os.unlink(path)

    def test_exon_interval_count(self):
        path = _write_temp(GFF3_BASIC)
        try:
            _, models = parse_gff(path, fmt='gff3')
            gm = models[0]
            self.assertEqual(len(gm.exon_intervals), 3)
        finally:
            os.unlink(path)

    def test_strand_plus(self):
        path = _write_temp(GFF3_BASIC)
        try:
            _, models = parse_gff(path, fmt='gff3')
            self.assertEqual(models[0].strand, '+')
        finally:
            os.unlink(path)

    def test_strand_minus(self):
        path = _write_temp(GFF3_MINUS)
        try:
            _, models = parse_gff(path, fmt='gff3')
            self.assertEqual(len(models), 1)
            self.assertEqual(models[0].strand, '-')
        finally:
            os.unlink(path)

    def test_partial_flag(self):
        path = _write_temp(GFF3_PARTIAL)
        try:
            _, models = parse_gff(path, fmt='gff3')
            gm = models[0]
            self.assertTrue(gm.is_partial_5)
        finally:
            os.unlink(path)

    def test_embedded_fasta(self):
        path = _write_temp(GFF3_WITH_FASTA)
        try:
            seqs, models = parse_gff(path, fmt='gff3')
            self.assertIn('TestSeq', seqs)
            self.assertGreater(len(seqs['TestSeq']), 0)
        finally:
            os.unlink(path)


# ─────────────────────────────────────────────────────────────────────────────
# GTF parsing tests
# ─────────────────────────────────────────────────────────────────────────────

class TestGTFParser(unittest.TestCase):

    def test_gtf_parse(self):
        path = _write_temp(GTF_BASIC, '.gtf')
        try:
            _, models = parse_gff(path, fmt='gtf')
            self.assertGreater(len(models), 0)
        finally:
            os.unlink(path)

    def test_gtf_gene_id(self):
        path = _write_temp(GTF_BASIC, '.gtf')
        try:
            _, models = parse_gff(path, fmt='gtf')
            gene_ids = [gm.gene_id for gm in models]
            self.assertIn('GTFGene1', gene_ids)
        finally:
            os.unlink(path)

    def test_gtf_cds_intervals(self):
        path = _write_temp(GTF_BASIC, '.gtf')
        try:
            _, models = parse_gff(path, fmt='gtf')
            gm = models[0]
            self.assertGreater(len(gm.cds_intervals), 0)
        finally:
            os.unlink(path)


# ─────────────────────────────────────────────────────────────────────────────
# GFF2 fallback tests
# ─────────────────────────────────────────────────────────────────────────────

class TestGFF2Parser(unittest.TestCase):

    def test_gff2_parse(self):
        path = _write_temp(GFF2_BASIC, '.gff')
        try:
            _, models = parse_gff(path, fmt='gff2')
            self.assertGreater(len(models), 0)
        finally:
            os.unlink(path)

    def test_gff2_gene_id(self):
        path = _write_temp(GFF2_BASIC, '.gff')
        try:
            _, models = parse_gff(path, fmt='gff2')
            self.assertEqual(models[0].gene_id, 'GFF2Gene1')
        finally:
            os.unlink(path)


# ─────────────────────────────────────────────────────────────────────────────
# Live data integration test
# ─────────────────────────────────────────────────────────────────────────────

class TestLiveGFF3Data(unittest.TestCase):

    GFF3_PATH = (
        '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
        'Augustus_ to Feature_file_to_genbank_file/Training_data/'
        'BAC_7C_17.gff3'
    )
    FASTA_PATH = (
        '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
        'Augustus_ to Feature_file_to_genbank_file/Training_data/'
        'BAC_7C17_complete_fasta.fasta'
    )

    @unittest.skipUnless(
        os.path.isfile(
            '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
            'Augustus_ to Feature_file_to_genbank_file/Training_data/'
            'BAC_7C_17.gff3'),
        'Training GFF3 not found'
    )
    def test_bac7c17_gff3_parses(self):
        seqs, models = parse_gff(self.GFF3_PATH,
                                  fasta_path=self.FASTA_PATH,
                                  fmt='gff3')
        self.assertGreater(len(models), 0)
        print(f'\n  [INFO] BAC_7C17 GFF3: {len(models)} gene model(s)')
        gene_ids = {gm.gene_id for gm in models}
        print(f'  [INFO] Gene IDs: {sorted(gene_ids)[:10]}')

    @unittest.skipUnless(
        os.path.isfile(
            '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
            'Augustus_ to Feature_file_to_genbank_file/Training_data/'
            'BAC_7C_17.gff3'),
        'Training GFF3 not found'
    )
    def test_hma4_genes_in_gff3(self):
        seqs, models = parse_gff(self.GFF3_PATH,
                                  fasta_path=self.FASTA_PATH,
                                  fmt='gff3')
        gene_ids = {gm.gene_id for gm in models}
        hma4_genes = {g for g in gene_ids if 'HMA4' in g}
        self.assertGreater(len(hma4_genes), 0,
                           f'Expected HMA4 genes, found: {gene_ids}')
        print(f'\n  [INFO] HMA4 genes in GFF3: {hma4_genes}')


if __name__ == '__main__':
    unittest.main(verbosity=2)
