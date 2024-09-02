from typing import Tuple, Dict, List, Optional
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt

# Known genes and their positions within the mitochondrial DNA sequence
KNOWN_GENES = {
    'MT-ND1': (3307, 4262),
    'MT-ND2': (4470, 5511),
    'MT-CO1': (5904, 7445),
    'MT-CO2': (7586, 8269),
    'MT-ATP6': (8527, 9207),
    'MT-ATP8': (8366, 8572),
    'MT-CYB': (14747, 15887),
    'MT-ND4': (10760, 12137),
    'MT-ND5': (12337, 14148)
}


def annotate_genes(sequence: SeqRecord, known_genes: Dict[str, Tuple[int, int]]) -> Dict[str, Tuple[str, int, int]]:
    """
    Annotate known genes in the sequence by extracting the subsequences
    that correspond to each gene's start and end positions.
    """
    annotations = {}
    for gene, (start, end) in known_genes.items():
        gene_seq = sequence.seq[start - 1:end]  # Extract the gene sequence
        annotations[gene] = (str(gene_seq), start, end)
    return annotations


def normalize_color(color: Tuple[int, int, int]) -> Tuple[float, float, float]:
    """
    Normalize an RGB tuple with values in the range 0-255 to a tuple with values in the range 0-1.
    """
    return tuple([c / 255 for c in color])


def plot_gene_annotation(sequence: SeqRecord, annotations: Dict[str, Tuple[str, int, int]], color: Tuple[int, int, int]) -> None:
    """
    Plot the annotated genes on the mitochondrial DNA sequence with the specified color.
    """
    color = normalize_color(color)

    fig, ax = plt.subplots(figsize=(14, 6))

    sequence_length = len(sequence)

    for gene, (seq, start, end) in annotations.items():
        ax.broken_barh([(start, end - start)], (10, 9), facecolors=[color], label=gene)
        ax.text((start + end) / 2, 14, gene, ha='center', va='center', rotation=90)

    ax.set_ylim(5, 25)
    ax.set_xlim(0, sequence_length)
    ax.set_yticks([15])
    ax.set_yticklabels(['Mitochondrial DNA'])
    ax.set_xlabel('Position in Sequence')
    ax.set_title(f'Annotated Mitochondrial Genes in {sequence.id}')
    ax.grid(True, linestyle='--', alpha=0.6)
    plt.show()


def gene_annotation_pipeline(sequences: List[SeqRecord], color: Tuple[int, int, int] = (0, 0, 255), start: Optional[int] = None, end: Optional[int] = None) -> None:
    """
    Gene annotation pipeline that processes a specified range of sequences, annotates known genes,
    and plots the results using the specified color.
    """
    if start is None and end is None:
        start = 0
        end = len(sequences)
    else:
        if start is None:
            start = 0
        if end is None:
            end = len(sequences)

    start = max(0, start)
    end = min(len(sequences), end)

    for sequence in sequences[start:end]:
        print(f"Processing sequence: {sequence.id}")
        annotations = annotate_genes(sequence, KNOWN_GENES)

        for gene, (seq, gene_start, gene_end) in annotations.items():
            print(f"{gene}: {seq[:30]}... (length: {len(seq)} bp)")

        plot_gene_annotation(sequence, annotations, color)
