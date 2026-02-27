"""
models.py
---------
Data classes used throughout the toolkit to represent gene models,
CDS intervals, and sequence metadata.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple


@dataclass
class CDSInterval:
    """
    A single CDS exon interval in 1-based genomic coordinates.

    Attributes
    ----------
    start : int
        Genomic start position (1-based, always <= stop, strand-independent).
    stop : int
        Genomic stop position (1-based, always >= start).
    phase : int
        GFF phase (0, 1, or 2). Number of bases at the start of this interval
        that belong to the previous codon. Used to derive codon_start for
        5'-partial CDS features.
    """
    start: int
    stop: int
    phase: int = 0


@dataclass
class ExonInterval:
    """
    A single mRNA exon interval (may include UTR).

    Attributes
    ----------
    start : int  Genomic start (1-based, <= stop).
    stop  : int  Genomic stop  (1-based, >= start).
    """
    start: int
    stop: int


@dataclass
class GeneModel:
    """
    Represents one predicted or annotated gene model with all information
    needed to write NCBI BankIt submission files.

    Coordinates
    -----------
    All coordinates are stored in 1-based, strand-independent form
    (start <= stop). The ``strand`` field carries directionality.
    The feature_writer converts these to BankIt convention (minus-strand
    features have first coord > second coord in the output table).

    Partiality
    ----------
    is_partial_5 : True when the 5' end of the CDS is incomplete
                   (no start codon). BankIt notation: ``<`` before the
                   5'-most coordinate of the CDS.
    is_partial_3 : True when the 3' end of the CDS is incomplete
                   (no stop codon). BankIt notation: ``>`` after the
                   3'-most coordinate of the CDS.
    """
    # Identifiers
    gene_id:       str
    transcript_id: str
    seq_id:        str            # matches FASTA header SeqID

    # Location
    strand:       str             # '+' or '-'
    gene_start:   int             # 1-based genomic start of gene span
    gene_stop:    int             # 1-based genomic stop  of gene span

    # Intervals (start <= stop, strand-independent)
    cds_intervals:  List[CDSInterval]   = field(default_factory=list)
    exon_intervals: List[ExonInterval]  = field(default_factory=list)  # mRNA exons

    # Partiality flags
    is_partial_5: bool = False
    is_partial_3: bool = False

    # Sequences (optional, may be filled by AUGUSTUS parser)
    coding_seq:  Optional[str] = None  # CDS nucleotide sequence (spliced)
    protein_seq: Optional[str] = None  # Protein sequence (no stop '*')

    # Score from AUGUSTUS prediction
    score: float = 0.0

    # Qualifiers for the BankIt feature table
    qualifiers: Dict[str, str] = field(default_factory=dict)
    # Common keys: 'gene', 'locus_tag', 'product', 'note', 'protein_id',
    #              'pseudo', 'codon_start'

    def cds_sorted_for_table(self) -> List[CDSInterval]:
        """
        Return CDS intervals sorted in transcript (5'→3') order for writing
        to the BankIt feature table.

        Plus strand  : ascending genomic coordinate (smallest first).
        Minus strand : descending genomic coordinate (largest first).
        """
        sorted_ivs = sorted(self.cds_intervals, key=lambda iv: iv.start)
        if self.strand == '-':
            sorted_ivs = sorted_ivs[::-1]
        return sorted_ivs

    def exon_sorted_for_table(self) -> List[ExonInterval]:
        """
        Same ordering logic for mRNA exon intervals.
        """
        sorted_ivs = sorted(self.exon_intervals, key=lambda iv: iv.start)
        if self.strand == '-':
            sorted_ivs = sorted_ivs[::-1]
        return sorted_ivs


@dataclass
class SequenceRecord:
    """
    A genomic sequence with optional BankIt-formatted metadata.

    Attributes
    ----------
    seq_id    : str  Identifier matching the FASTA header (no spaces).
    sequence  : str  Full nucleotide sequence (uppercase or lowercase).
    organism  : str  Scientific name, e.g. "Arabidopsis halleri subsp. halleri".
    mol_type  : str  Molecule type, typically "genomic DNA".
    modifiers : dict Additional BankIt source modifiers:
                     isolate, strain, clone, country, collection_date, etc.
    description : str Free-text description for the FASTA defline.
    """
    seq_id:      str
    sequence:    str
    organism:    str  = "Organism name not specified"
    mol_type:    str  = "genomic DNA"
    modifiers:   Dict[str, str] = field(default_factory=dict)
    description: str = "genomic sequence"
