# AUGUSTUS BankIt Toolkit

Convert **AUGUSTUS gene predictions** and existing **GenBank / GFF / GTF** annotations into **NCBI BankIt-compatible** submission files — fully following NCBI rules and the INSDC feature table specification.

---

## Overview

Submitting genomic sequences with gene annotations to NCBI GenBank via BankIt requires three files:

| File | Description |
|------|-------------|
| `*_nucleotide.fasta` | Genomic FASTA with BankIt source-modifier deflines |
| `*_feature_table.tbl` | 5-column tab-delimited feature table (INSDC format) |
| `*_protein.fasta` | Protein sequences (reference/validation) |

This toolkit provides **two conversion tools** and one **validator**:

```
augustus2bankit.py   →  AUGUSTUS HTML/GFF output  →  BankIt files
genbank2bankit.py    →  GenBank / GFF3 / GTF      →  BankIt files
validate_bankit.py   →  Check output files before submission
```

---

## Tools

### 1. `augustus2bankit.py` — AUGUSTUS → BankIt

Parses the HTML output from the [AUGUSTUS web server](https://bioinf.uni-greifswald.de/augustus/)
(or plain AUGUSTUS GFF text) and generates BankIt submission files.

**What it extracts from AUGUSTUS HTML:**
- Gene coordinates (all exons, strand)
- CDS intervals and phases
- Coding nucleotide sequences (inline in HTML)
- Protein sequences (inline in HTML)
- Partiality flags (no start/stop codon → `<`/`>` notation)

**What you provide:**
- The original genomic FASTA that was submitted to AUGUSTUS
- Organism name and BankIt source modifiers

```bash
python augustus2bankit.py \
    --augustus  Augustus_BAC_7C17.html \
    --fasta     BAC_7C17.fasta \
    --organism  "Arabidopsis halleri subsp. halleri" \
    --isolate   "BAC_7C17" \
    --product   "Zn/Cd P(IB)-type ATPase" \
    --outdir    bankit_output/
```

---

### 2. `genbank2bankit.py` — GenBank/GFF/GTF → BankIt  *(reverse tool)*

Extracts gene models from existing annotation files and reformats them
as BankIt submission files. Supports:

| Format | File extensions | Notes |
|--------|----------------|-------|
| GenBank flat file | `.gb`, `.genbank`, `.gbk`, `.txt` | Full sequence + features |
| GFF3 | `.gff3`, `.gff` | Embedded or separate FASTA |
| GTF  | `.gtf` | Separate FASTA required |
| GFF2 | `.gff`, `.gff2` | Fallback parsing |

```bash
# From GenBank file (extract only HMA4 genes)
python genbank2bankit.py \
    --input     BAC_7C17_training_data_from_NCBI.gb \
    --format    genbank \
    --organism  "Arabidopsis halleri subsp. halleri" \
    --clone     "BAC_7C17" \
    --gene-list "AhHMA4-1,AhHMA4-2" \
    --outdir    bankit_output/

# From GFF3 + separate FASTA
python genbank2bankit.py \
    --input     BAC_7C_17.gff3 \
    --fasta     BAC_7C17_complete_fasta.fasta \
    --organism  "Arabidopsis halleri subsp. halleri" \
    --outdir    bankit_output/
```

---

### 3. `validate_bankit.py` — Pre-submission Validator

Runs a comprehensive validation suite before you upload to BankIt.

```bash
python validate_bankit.py \
    --tbl   output_feature_table.tbl \
    --nuc   output_nucleotide.fasta \
    --prot  output_protein.fasta \
    --report validation_report.txt
```

**Checks performed:**

| Category | Checks |
|----------|--------|
| Feature table | Tab structure, INSDC feature keys, coordinate range, required qualifiers |
| Nucleotide FASTA | `[organism=]` and `[mol_type=]` present, IUPAC characters, no duplicates |
| Protein FASTA | Starts with M, no internal stop codons, valid amino acids |
| Cross-file | SeqID matching, CDS length divisible by 3, translation consistency |

**Exit codes:** `0` = pass, `1` = errors found, `2` = file not found.

---

## Installation

```bash
# Clone the repository
git clone https://github.com/ap1438/annotation2bankit.git
cd annotation2bankit

# No external dependencies required (Python 3.6+ standard library only)
# Optional: install for test suite
pip install pytest

# Or install as a package
pip install .
```

---

## Quick Start

### Example 1: AUGUSTUS web report → BankIt (HMA4 genes, unknown ecotype)

```bash
# 1. Run AUGUSTUS on your sequences at https://bioinf.uni-greifswald.de/augustus/
#    Save the HTML result page (File → Save As...).

# 2. Convert
python augustus2bankit.py \
    --augustus  "Augustus_BAC_7C17.html" \
    --fasta     "BAC_7C17.fasta" \
    --organism  "Arabidopsis halleri subsp. halleri" \
    --isolate   "BAC_7C17" \
    --sub-species "halleri" \
    --country   "Germany: Harz mountains, Langelsheim" \
    --date      "2007" \
    --product   "Zn/Cd P(IB)-type ATPase" \
    --locus-tag "AhHMA4" \
    --outdir    "bankit_BAC_7C17/" \
    --prefix    "BAC_7C17_hma4" \
    --primary-only

# 3. Validate
python validate_bankit.py \
    --tbl   bankit_BAC_7C17/BAC_7C17_hma4_feature_table.tbl \
    --nuc   bankit_BAC_7C17/BAC_7C17_hma4_nucleotide.fasta \
    --prot  bankit_BAC_7C17/BAC_7C17_hma4_protein.fasta

# 4. Upload the .tbl and _nucleotide.fasta to https://www.ncbi.nlm.nih.gov/WebSub/
```

### Example 2: GenBank file → extract HMA4 → new BankIt files

```bash
python genbank2bankit.py \
    --input     EU382073.gb \
    --organism  "Arabidopsis halleri subsp. halleri" \
    --gene-list "AhHMA4-1,AhHMA4-2" \
    --clone     "BAC_7C17" \
    --outdir    bankit_hma4/ \
    --prefix    hma4_resubmission
```

---

## NCBI BankIt Rules Implemented

This toolkit strictly follows the rules documented in:
- [BankIt Feature Table Format](https://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html)
- [INSDC Feature Table Documentation](https://www.insdc.org/submitting-standards/feature-table-documentation/)
- [BankIt Nucleotide FASTA Submission Help](https://www.ncbi.nlm.nih.gov/WebSub/html/help/nucleotide-submissions.html)

### Feature Table Format Rules

```
>Feature SeqID
start   stop    gene
                gene    GeneName
start   stop    mRNA
cont1   cont1
cont2   cont2
                product protein product name
start   stop    CDS
cont1   cont1
cont2   cont2
                product protein product name
                gene    GeneName
```

Key rules:
- **Tab-delimited**: start `\t` stop `\t` feature_key
- **Qualifiers**: exactly **3 leading tabs** then qualifier_key `\t` qualifier_value
- **Continuation intervals**: just start `\t` stop (no feature key)
- **Minus strand**: first coordinate > second coordinate
- **Multi-exon minus strand**: listed in 5'→3' transcript order (descending genomic position)
- **5'-partial**: `<` before the 5'-most coordinate
- **3'-partial**: `>` before the 3'-most coordinate
- **SeqID**: must exactly match the FASTA defline identifier (no spaces)

### Nucleotide FASTA Format Rules

```
>SeqID [organism=Scientific name] [mol_type=genomic DNA] [isolate=name] description
ACGTACGTACGT...
ACGTACGTACGT...
```

Required modifiers: `[organism=]` `[mol_type=]`
Sequence lines: 60 characters (NCBI standard).

### Supported Feature Qualifiers

| Qualifier | Where | Notes |
|-----------|-------|-------|
| `gene` | gene, mRNA, CDS | Gene symbol |
| `locus_tag` | gene, mRNA, CDS | Unique locus identifier |
| `product` | mRNA, CDS | Protein product description (**required for CDS**) |
| `protein_id` | CDS | Format: `accession.version` |
| `codon_start` | CDS | 1, 2, or 3 (for 5'-partial CDS) |
| `note` | any | Free text note |
| `pseudo` | gene | Boolean – marks pseudogene |

---

## Project Structure

```
augustus_bankit_toolkit/
├── README.md                  ← This file
├── LICENSE
├── requirements.txt
├── setup.py
│
├── augustus2bankit.py         ← Tool 1: AUGUSTUS → BankIt
├── genbank2bankit.py          ← Tool 2: GenBank/GFF/GTF → BankIt
├── validate_bankit.py         ← Tool 3: Validation
│
├── lib/
│   ├── __init__.py
│   ├── models.py              ← GeneModel, CDSInterval, SequenceRecord
│   ├── fasta_utils.py         ← FASTA I/O, translation, reverse complement
│   ├── augustus_parser.py     ← AUGUSTUS HTML/GFF parser
│   ├── genbank_parser.py      ← GenBank flat file parser
│   ├── gff_parser.py          ← GFF2/GFF3/GTF parser
│   ├── feature_writer.py      ← BankIt feature table + FASTA writer
│   └── ncbi_validator.py      ← Validation suite

```

---

## Command-Line Reference

### `augustus2bankit.py`

```
Required:
  --augustus FILE     AUGUSTUS HTML or plain GFF output
  --fasta    FILE     Genomic FASTA (submitted to AUGUSTUS)
  --organism STR      Scientific name

Source modifiers:
  --mol-type    STR   Molecule type (default: "genomic DNA")
  --isolate     STR   Isolate name
  --strain      STR   Strain name
  --clone       STR   Clone identifier
  --country     STR   "Country:Region" format
  --date        STR   Collection date (e.g. "2007")
  --sub-species STR   Sub-species
  --genotype    STR   Genotype

Annotation overrides:
  --gene-name  STR    Gene name for all models
  --product    STR    Product description for all CDS
  --locus-tag  STR    Locus tag prefix (models get -1, -2, ...)
  --gene-list  FILE   List of gene IDs to include

Output:
  --outdir  DIR       Output directory (default: .)
  --prefix  STR       File prefix (default: output)

Flags:
  --primary-only      Keep only first isoform per gene
  --no-validate       Skip validation
  --verbose           Verbose output
```

### `genbank2bankit.py`

```
Required:
  --input    FILE     Annotation file (GenBank/GFF3/GTF/GFF2)
  --organism STR      Scientific name

Input options:
  --fasta    FILE     Separate FASTA (for GFF without embedded sequences)
  --format   STR      auto|genbank|gff3|gtf|gff2 (default: auto)

Filtering:
  --gene-list STR     Comma-separated gene names or file path

(Same source modifiers and output options as augustus2bankit.py)
```

### `validate_bankit.py`

```
Required:
  --tbl  FILE    Feature table (.tbl)
  --nuc  FILE    Nucleotide FASTA

Optional:
  --prot   FILE  Protein FASTA
  --report FILE  Save report to file
  --verbose      Detailed output

Exit codes: 0=pass, 1=errors, 2=file missing
```

---

## Troubleshooting

### "No gene models found in AUGUSTUS output"

- Make sure the HTML file was saved completely (not just the visible portion)
- Try saving the AUGUSTUS result page as "Web Page, Complete" or "HTML only"
- The HTML must contain a `<pre class="result">` block

### "SeqID in feature table has no matching sequence in FASTA"

- The sequence names in your FASTA must exactly match those in the AUGUSTUS output
- Check the first token of each FASTA defline (before the first space)

### "CDS total length X bp is not divisible by 3"

- This is a genuine annotation error — coordinate boundaries may be off by one
- Check your AUGUSTUS settings (species model, min exon prob)
- For partial genes, ensure `--partial` handling is correct in AUGUSTUS

### Feature table shows wrong strand

- AUGUSTUS minus-strand genes have CDS on the reverse complement
- The BankIt feature table uses `c1 > c2` convention for minus strand
- Verify AUGUSTUS GFF shows `-` in column 7

---

## Citation

If you use this toolkit, please cite:

> AUGUSTUS: Stanke, M., Diekhans, M., Baertsch, R. & Haussler, D. (2008).
> Using native and syntenically mapped cDNA alignments to improve de novo gene finding.
> *Bioinformatics* 24, 637–644. doi:10.1093/bioinformatics/btn013

For the HMA4 gene annotation (training data):

> Hanikenne, M. et al. (2008). Evolution of metal hyperaccumulation required
> cis-regulatory changes and triplication of HMA4.
> *Nature* 453, 391–395. doi:10.1038/nature06877

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contributing

Pull requests are welcome. Please:
1. Add tests for any new features
2. Run `bash tests/run_tests.sh` before submitting
3. Follow PEP 8 code style
4. Update this README for any new CLI options
