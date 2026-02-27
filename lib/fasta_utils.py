"""
fasta_utils.py
--------------
FASTA file reading/writing, sequence manipulation, and CDS translation.

NCBI BankIt FASTA rules implemented here:
  - Defline: ">SeqID [modifier=value] ... description"
  - SeqID   : no spaces (first whitespace-delimited token after '>')
  - Sequence : wrapped at 60 characters per line (NCBI default)
  - Modifiers: bracketed key=value pairs on the defline
  - Required modifiers for genomic submissions:
        [organism=...] [mol_type=genomic DNA]
  - Common optional modifiers:
        [isolate=] [strain=] [clone=] [country=] [collection_date=]
        [sub_species=] [genotype=]
"""

import re
from typing import Dict, Optional, List, Tuple

# ─────────────────────────────────────────────────────────────────────────────
# Genetic code table 1 (standard / universal)
# ─────────────────────────────────────────────────────────────────────────────
CODON_TABLE: Dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

RC_TABLE = str.maketrans('ACGTacgt', 'TGCAtgca')


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(RC_TABLE)[::-1]


def translate_cds(cds_seq: str, codon_start: int = 1) -> str:
    """
    Translate a CDS nucleotide sequence to a protein sequence.

    Parameters
    ----------
    cds_seq     : str  Nucleotide sequence (5'→3' of the coding strand).
    codon_start : int  1, 2, or 3 (NCBI codon_start qualifier).
                       1 = start at position 1 (normal).
                       2 = skip first base (5' partial, reading frame +1).
                       3 = skip first two bases (5' partial, reading frame +2).

    Returns
    -------
    str : Protein sequence WITHOUT trailing stop codon ('*').
    """
    seq = cds_seq.upper().replace('\n', '').replace(' ', '')
    offset = codon_start - 1  # codon_start=1 → offset=0
    seq = seq[offset:]
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


def read_fasta(filepath: str) -> Dict[str, str]:
    """
    Parse a FASTA file and return a dict mapping SeqID → sequence.

    The SeqID is the first whitespace-delimited token after '>'.
    The sequence is returned as an uppercase string with no whitespace.
    """
    sequences: Dict[str, str] = {}
    current_id: Optional[str] = None
    current_seq: List[str] = []

    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq).upper()
                # SeqID = first token after '>'
                header = line[1:]
                current_id = header.split()[0] if header.split() else header
                current_seq = []
            elif current_id is not None:
                current_seq.append(line.strip())

    if current_id is not None:
        sequences[current_id] = ''.join(current_seq).upper()

    return sequences


def format_bankit_defline(
    seq_id: str,
    organism: str,
    mol_type: str = 'genomic DNA',
    description: str = 'genomic sequence',
    **modifiers: str
) -> str:
    """
    Build a BankIt-compatible FASTA defline.

    NCBI BankIt format:
        >SeqID [organism=value] [mol_type=value] [key=value ...] description

    Required modifiers (always written first):
        organism, mol_type

    Common optional modifiers (written in this order when present):
        isolate, strain, clone, sub_species, country, collection_date,
        genotype, cultivar, ecotype, sex, dev_stage

    Parameters
    ----------
    seq_id      : Identifier (no spaces).
    organism    : Scientific name.
    mol_type    : Molecule type (default: 'genomic DNA').
    description : Free-text description appended after modifiers.
    **modifiers : Additional key=value modifier pairs.

    Returns
    -------
    str : Full defline string (including the leading '>').
    """
    # Ordered optional modifier keys
    ordered_keys = [
        'isolate', 'strain', 'clone', 'sub_species', 'country',
        'collection_date', 'genotype', 'cultivar', 'ecotype',
        'sex', 'dev_stage', 'note',
    ]
    parts = [f'>{seq_id}']
    parts.append(f'[organism={organism}]')
    parts.append(f'[mol_type={mol_type}]')
    for key in ordered_keys:
        if key in modifiers and modifiers[key]:
            parts.append(f'[{key}={modifiers[key]}]')
    # Any remaining modifiers not in the ordered list
    for key, val in modifiers.items():
        if key not in ordered_keys and val:
            parts.append(f'[{key}={val}]')
    if description:
        parts.append(description)
    return ' '.join(parts)


def write_fasta(
    sequences: Dict[str, str],
    filepath: str,
    width: int = 60,
    deflines: Optional[Dict[str, str]] = None,
) -> None:
    """
    Write sequences to a FASTA file.

    Parameters
    ----------
    sequences : dict  {seq_id: sequence}
    filepath  : str   Output file path.
    width     : int   Characters per sequence line (default 60, NCBI standard).
    deflines  : dict  Optional {seq_id: full_defline_without_arrow}.
                      If None, bare '>seq_id' is used.
    """
    with open(filepath, 'w') as fh:
        for seq_id, seq in sequences.items():
            if deflines and seq_id in deflines:
                fh.write(f'>{deflines[seq_id]}\n')
            else:
                fh.write(f'>{seq_id}\n')
            seq_upper = seq.upper()
            for i in range(0, len(seq_upper), width):
                fh.write(seq_upper[i:i + width] + '\n')


def extract_cds_from_genome(
    genomic_seq: str,
    intervals: list,   # list of (start, stop) tuples, 1-based, start<=stop
    strand: str,
) -> str:
    """
    Extract and splice a CDS from a genomic sequence.

    Parameters
    ----------
    genomic_seq : str  Full genomic sequence (0-indexed internally).
    intervals   : list List of (start, stop) tuples, 1-based coordinates,
                       start <= stop regardless of strand.
    strand      : str  '+' or '-'.

    Returns
    -------
    str : Spliced CDS in 5'→3' coding-strand orientation (uppercase).

    Notes
    -----
    Intervals should be sorted in ascending genomic order before calling.
    This function will sort them ascending and then reverse-complement the
    entire spliced sequence for minus-strand genes.
    """
    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    spliced = []
    for start, stop in sorted_ivs:
        # Convert to 0-based Python slice
        spliced.append(genomic_seq[start - 1: stop].upper())
    cds = ''.join(spliced)
    if strand == '-':
        cds = reverse_complement(cds)
    return cds


def validate_fasta_file(filepath: str) -> Tuple[bool, List[str]]:
    """
    Basic FASTA format validation for BankIt submissions.

    Checks:
    1. File is not empty.
    2. Every record starts with '>'.
    3. SeqID has no spaces.
    4. Sequence contains only IUPAC nucleotide characters.
    5. No duplicate SeqIDs.

    Returns
    -------
    (is_valid: bool, messages: list of str)
    """
    errors: List[str] = []
    seen_ids: set = set()
    current_id: Optional[str] = None
    has_sequence = False
    iupac = set('ACGTURYSWKMBDHVNacgturyswkmbdhvn')

    try:
        with open(filepath, 'r') as fh:
            lines = fh.readlines()
    except FileNotFoundError:
        return False, [f'File not found: {filepath}']

    if not lines:
        return False, ['FASTA file is empty.']

    for lineno, line in enumerate(lines, 1):
        line = line.rstrip('\n')
        if not line:
            continue
        if line.startswith('>'):
            tokens = line[1:].split()
            if not tokens:
                errors.append(f'Line {lineno}: Empty header after ">".')
                current_id = None
                continue
            current_id = tokens[0]
            if ' ' in current_id:
                errors.append(
                    f'Line {lineno}: SeqID "{current_id}" contains a space.')
            if current_id in seen_ids:
                errors.append(
                    f'Line {lineno}: Duplicate SeqID "{current_id}".')
            seen_ids.add(current_id)
        else:
            has_sequence = True
            invalid_chars = set(line.strip()) - iupac
            if invalid_chars:
                errors.append(
                    f'Line {lineno}: Invalid nucleotide characters: '
                    f'{invalid_chars}')

    if not has_sequence:
        errors.append('No sequence data found in FASTA file.')

    return (len(errors) == 0), errors
