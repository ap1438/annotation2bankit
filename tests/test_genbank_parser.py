#!/usr/bin/env python3
"""
test_genbank_parser.py
----------------------
Tests for the GenBank flat file parser.

Tests verify:
  1. LOCUS, DEFINITION, ACCESSION, ORGANISM parsing
  2. CDS feature extraction (location, qualifiers, translation)
  3. mRNA feature extraction
  4. Complement strand handling
  5. Multi-exon join() locations
  6. Partial features (<, >)
  7. Multi-record GenBank files
  8. Integration: parsing training data file

Run:
    python -m pytest tests/test_genbank_parser.py -v
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lib.genbank_parser import parse_genbank, _parse_location
from lib.models import GeneModel


# ─────────────────────────────────────────────────────────────────────────────
# Minimal GenBank records for testing
# ─────────────────────────────────────────────────────────────────────────────

SIMPLE_GB = """\
LOCUS       TEST001               1000 bp    DNA     linear   PLN 01-JAN-2024
DEFINITION  Test sequence, genomic sequence.
ACCESSION   TEST001
VERSION     TEST001.1
KEYWORDS    .
SOURCE      Test organism
  ORGANISM  Test organism
            Eukaryota.
FEATURES             Location/Qualifiers
     source          1..1000
                     /organism="Test organism"
                     /mol_type="genomic DNA"
     gene            100..900
                     /gene="TestGene1"
     mRNA            join(100..300,500..700,800..900)
                     /gene="TestGene1"
                     /product="test protein"
     CDS             join(100..300,500..700,800..900)
                     /gene="TestGene1"
                     /product="test protein"
                     /codon_start=1
                     /translation="MAAAA"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
       61 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
      121 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""

COMPLEMENT_GB = """\
LOCUS       TEST002               1000 bp    DNA     linear   PLN 01-JAN-2024
DEFINITION  Test complement sequence.
ACCESSION   TEST002
VERSION     TEST002.1
KEYWORDS    .
SOURCE      Test organism
  ORGANISM  Test organism
            Eukaryota.
FEATURES             Location/Qualifiers
     source          1..1000
                     /organism="Test organism"
                     /mol_type="genomic DNA"
     gene            complement(100..900)
                     /gene="TestGene2"
     CDS             complement(join(100..300,500..700,800..900))
                     /gene="TestGene2"
                     /product="complement protein"
                     /translation="MBBBB"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""

PARTIAL_GB = """\
LOCUS       TEST003               1000 bp    DNA     linear   PLN 01-JAN-2024
DEFINITION  Partial test sequence.
ACCESSION   TEST003
VERSION     TEST003.1
KEYWORDS    .
SOURCE      Test organism
  ORGANISM  Test organism
            Eukaryota.
FEATURES             Location/Qualifiers
     source          1..1000
                     /organism="Test organism"
                     /mol_type="genomic DNA"
     gene            <1..>500
                     /locus_tag="TL001"
     CDS             <1..>500
                     /locus_tag="TL001"
                     /product="partial protein"
                     /codon_start=2
                     /translation="RKKK"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""

MULTI_RECORD_GB = SIMPLE_GB + COMPLEMENT_GB


def _write_temp_gb(content: str) -> str:
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gb', delete=False) as fh:
        fh.write(content)
        return fh.name


# ─────────────────────────────────────────────────────────────────────────────
# Location parser unit tests
# ─────────────────────────────────────────────────────────────────────────────

class TestLocationParser(unittest.TestCase):

    def test_simple_range(self):
        ivs, strand, p5, p3 = _parse_location('100..500')
        self.assertEqual(ivs, [(100, 500)])
        self.assertEqual(strand, '+')
        self.assertFalse(p5)
        self.assertFalse(p3)

    def test_complement(self):
        ivs, strand, p5, p3 = _parse_location('complement(100..500)')
        self.assertEqual(ivs, [(100, 500)])
        self.assertEqual(strand, '-')

    def test_join(self):
        ivs, strand, p5, p3 = _parse_location('join(100..200,300..500)')
        self.assertEqual(ivs, [(100, 200), (300, 500)])
        self.assertEqual(strand, '+')

    def test_complement_join(self):
        ivs, strand, p5, p3 = _parse_location(
            'complement(join(100..200,300..500))')
        self.assertEqual(ivs, [(100, 200), (300, 500)])
        self.assertEqual(strand, '-')

    def test_partial_5(self):
        ivs, strand, p5, p3 = _parse_location('<1..500')
        self.assertTrue(p5)
        self.assertFalse(p3)

    def test_partial_3(self):
        ivs, strand, p5, p3 = _parse_location('100..>500')
        self.assertFalse(p5)
        self.assertTrue(p3)

    def test_partial_both(self):
        ivs, strand, p5, p3 = _parse_location('<1..>500')
        self.assertTrue(p5)
        self.assertTrue(p3)


# ─────────────────────────────────────────────────────────────────────────────
# GenBank parser tests
# ─────────────────────────────────────────────────────────────────────────────

class TestGenBankParserBasic(unittest.TestCase):

    def test_simple_record_parsed(self):
        path = _write_temp_gb(SIMPLE_GB)
        try:
            results = parse_genbank(path)
            self.assertEqual(len(results), 1)
            rec, models = results[0]
            self.assertEqual(rec.seq_id, 'TEST001')
        finally:
            os.unlink(path)

    def test_sequence_extracted(self):
        path = _write_temp_gb(SIMPLE_GB)
        try:
            results = parse_genbank(path)
            rec, _ = results[0]
            self.assertGreater(len(rec.sequence), 0)
        finally:
            os.unlink(path)

    def test_organism_extracted(self):
        path = _write_temp_gb(SIMPLE_GB)
        try:
            results = parse_genbank(path)
            rec, _ = results[0]
            self.assertIn('Test', rec.organism)
        finally:
            os.unlink(path)

    def test_gene_model_extracted(self):
        path = _write_temp_gb(SIMPLE_GB)
        try:
            results = parse_genbank(path)
            _, models = results[0]
            self.assertGreater(len(models), 0)
            gm = models[0]
            self.assertIn('TestGene', gm.gene_id)
        finally:
            os.unlink(path)

    def test_cds_intervals(self):
        path = _write_temp_gb(SIMPLE_GB)
        try:
            results = parse_genbank(path)
            _, models = results[0]
            gm = models[0]
            self.assertEqual(len(gm.cds_intervals), 3)
        finally:
            os.unlink(path)

    def test_product_qualifier(self):
        path = _write_temp_gb(SIMPLE_GB)
        try:
            results = parse_genbank(path)
            _, models = results[0]
            gm = models[0]
            self.assertEqual(gm.qualifiers.get('product'), 'test protein')
        finally:
            os.unlink(path)

    def test_translation_stored(self):
        path = _write_temp_gb(SIMPLE_GB)
        try:
            results = parse_genbank(path)
            _, models = results[0]
            gm = models[0]
            self.assertIsNotNone(gm.protein_seq)
        finally:
            os.unlink(path)


class TestGenBankParserStrand(unittest.TestCase):

    def test_complement_strand(self):
        path = _write_temp_gb(COMPLEMENT_GB)
        try:
            results = parse_genbank(path)
            _, models = results[0]
            self.assertGreater(len(models), 0)
            gm = models[0]
            self.assertEqual(gm.strand, '-')
        finally:
            os.unlink(path)

    def test_complement_intervals(self):
        path = _write_temp_gb(COMPLEMENT_GB)
        try:
            results = parse_genbank(path)
            _, models = results[0]
            gm = models[0]
            self.assertEqual(len(gm.cds_intervals), 3)
            # All intervals should have start <= stop (genomic convention)
            for iv in gm.cds_intervals:
                self.assertLessEqual(iv.start, iv.stop)
        finally:
            os.unlink(path)


class TestGenBankParserPartial(unittest.TestCase):

    def test_partial_flags(self):
        path = _write_temp_gb(PARTIAL_GB)
        try:
            results = parse_genbank(path)
            _, models = results[0]
            self.assertGreater(len(models), 0)
            gm = models[0]
            self.assertTrue(gm.is_partial_5)
            self.assertTrue(gm.is_partial_3)
        finally:
            os.unlink(path)

    def test_locus_tag_extracted(self):
        path = _write_temp_gb(PARTIAL_GB)
        try:
            results = parse_genbank(path)
            _, models = results[0]
            gm = models[0]
            self.assertEqual(gm.qualifiers.get('locus_tag'), 'TL001')
        finally:
            os.unlink(path)


class TestGenBankParserMultiRecord(unittest.TestCase):

    def test_multi_record_count(self):
        path = _write_temp_gb(MULTI_RECORD_GB)
        try:
            results = parse_genbank(path)
            self.assertEqual(len(results), 2)
        finally:
            os.unlink(path)

    def test_multi_record_seq_ids(self):
        path = _write_temp_gb(MULTI_RECORD_GB)
        try:
            results = parse_genbank(path)
            ids = [rec.seq_id for rec, _ in results]
            self.assertIn('TEST001', ids)
            self.assertIn('TEST002', ids)
        finally:
            os.unlink(path)


class TestLiveGenBankData(unittest.TestCase):
    """Integration test against real training data."""

    GB_PATH = (
        '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
        'Augustus_ to Feature_file_to_genbank_file/Training_data/'
        'BAC_7C17_training_data_from_NCBI.gb'
    )

    @unittest.skipUnless(
        os.path.isfile(
            '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
            'Augustus_ to Feature_file_to_genbank_file/Training_data/'
            'BAC_7C17_training_data_from_NCBI.gb'),
        'Training data not found'
    )
    def test_bac7c17_gb_parses(self):
        results = parse_genbank(self.GB_PATH)
        self.assertGreater(len(results), 0)
        rec, models = results[0]
        print(f'\n  [INFO] GenBank: SeqID={rec.seq_id}, '
              f'len={len(rec.sequence)} bp, '
              f'{len(models)} gene model(s)')
        self.assertEqual(rec.seq_id, 'EU382073')
        self.assertGreater(len(models), 0)

    @unittest.skipUnless(
        os.path.isfile(
            '/prj/pflaphy-pacbio/Claud_AI/HMA4/'
            'Augustus_ to Feature_file_to_genbank_file/Training_data/'
            'BAC_7C17_training_data_from_NCBI.gb'),
        'Training data not found'
    )
    def test_hma4_genes_found(self):
        results = parse_genbank(self.GB_PATH)
        _, models = results[0]
        gene_ids = {gm.gene_id for gm in models}
        # AhHMA4-1 and AhHMA4-2 should be present
        hma4 = {g for g in gene_ids if 'HMA4' in g or 'hma4' in g.lower()}
        print(f'\n  [INFO] HMA4 genes found: {hma4}')
        self.assertGreater(len(hma4), 0, 'Expected at least one HMA4 gene')


if __name__ == '__main__':
    unittest.main(verbosity=2)
