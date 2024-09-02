import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
from Bio import SeqIO
from typing import Optional, Tuple


def preprocess_sequences(fasta_file: str, output_file: str) -> None:
    """
    Preprocess sequences to ensure they are of the same length.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))

    min_length = min(len(record.seq) for record in records)

    cropped_records = []
    for record in records:
        cropped_seq = record.seq[:min_length]
        cropped_record = SeqRecord(cropped_seq, id=record.id,
                                   description=record.description)
        cropped_records.append(cropped_record)

    with open(output_file, "w") as output_handle:
        SeqIO.write(cropped_records, output_handle, "fasta")


def load_alignment(fasta_file: str) -> MultipleSeqAlignment:
    """
    Load a multiple sequence alignment from a FASTA file.
    """
    temp_fasta_file = "preprocessed_conserved" + fasta_file
    preprocess_sequences(fasta_file, temp_fasta_file)

    alignment = AlignIO.read(temp_fasta_file, "fasta")
    return alignment


def calculate_conservation(alignment: MultipleSeqAlignment) -> np.ndarray:
    """
    Calculate conservation at each position in the alignment.
    """
    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()

    conservation_scores = np.zeros(alignment_length)

    for i in range(alignment_length):
        column = alignment[:, i]
        most_common_residue = max(set(column), key=column.count)
        conservation_scores[i] = column.count(most_common_residue) / num_sequences

    return conservation_scores


def plot_conserved_regions(conservation_scores: np.ndarray, color: Tuple[int, int, int]) -> None:
    """
    Plot the conservation scores across the sequence alignment.
    """
    normalized_color = tuple([c / 255 for c in color])

    plt.figure(figsize=(14, 6))
    sns.lineplot(x=range(len(conservation_scores)), y=conservation_scores, marker="o", color=normalized_color)
    plt.title("Conserved Regions in Sequence Alignment")
    plt.xlabel("Position in Alignment")
    plt.ylabel("Conservation Score")
    plt.ylim(0, 1)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()


def highlight_conserved_regions(conservation_scores: np.ndarray, threshold: float, color: Tuple[int, int, int]) -> None:
    """
    Highlight regions with conservation scores above the threshold.
    """
    normalized_color = tuple([c / 255 for c in color])

    plt.figure(figsize=(14, 6))
    x_values = range(len(conservation_scores))

    sns.lineplot(x=x_values, y=conservation_scores, marker="o", color=normalized_color, label="Conservation Score")

    plt.fill_between(x_values, conservation_scores, where=(conservation_scores >= threshold),
                     color=normalized_color, alpha=0.4, label=f"Conserved Regions (Score >= {threshold})")

    plt.title(f"Conserved Regions (Threshold: {threshold})")
    plt.xlabel("Position in Alignment")
    plt.ylabel("Conservation Score")
    plt.ylim(0, 1)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()

    plt.show()


def analyze_conserved_regions(fasta_file: str, threshold: float = 0.9,
                              color: Tuple[int, int, int] = (0, 0, 255),
                              start: Optional[int] = None, end: Optional[int]
                              = None) -> None:
    """
    Analyze and plot conserved regions from a sequence alignment.

    Parameters:
    - fasta_file: Path to the input FASTA file.
    - threshold: A float representing the conservation score threshold for highlighting regions (default is 0.9).
    - color: A tuple representing an RGB color value for the plot and highlighted regions (default is blue).
    - start: The start index of the sequence range to be processed (optional).
    - end: The end index of the sequence range to be processed (optional).
    """
    alignment = load_alignment(fasta_file)

    if start is None:
        start = 0
    if end is None:
        end = len(alignment)

    alignment = alignment[start:end]

    conservation_scores = calculate_conservation(alignment)

    plot_conserved_regions(conservation_scores, color)

    highlight_conserved_regions(conservation_scores, threshold, color)
