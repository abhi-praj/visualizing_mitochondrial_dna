import os

import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
from Bio import SeqIO
from typing import Tuple


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
    temp_fasta_file = "preprocessed_conserved_" + os.path.basename(fasta_file)
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
        conservation_scores[i] = column.count(
            most_common_residue) / num_sequences

    return conservation_scores


def plot_conserved_regions(conservation_scores: np.ndarray,
                           color: Tuple[int, int, int],
                           output_path: str) -> None:
    """
    Plot the conservation scores across the sequence alignment and save the plot as an image.
    """
    normalized_color = tuple(c / 255 for c in color)

    plt.figure(figsize=(14, 6))
    sns.lineplot(x=range(len(conservation_scores)), y=conservation_scores,
                 marker="o", color=normalized_color)
    plt.title("Conserved Regions in Sequence Alignment")
    plt.xlabel("Position in Alignment")
    plt.ylabel("Conservation Score")
    plt.ylim(0, 1)
    plt.grid(True, linestyle='--', alpha=0.6)

    plt.savefig(output_path)
    plt.close()


def analyze_conserved_regions(fasta_file: str, output_path_conserved: str,
                              threshold: float = 0.9,
                              color: Tuple[int, int, int] = (
                              0, 0, 0)) -> None:
    """
    Analyze and plot conserved regions from a sequence alignment, and save the plot as an image.

    Parameters:
    - fasta_file: Path to the input FASTA file.
    - output_path_conserved: Path to save the conserved regions plot.
    - threshold: Conservation score threshold for highlighting regions (default is 0.9).
    - color: RGB color value for the plot (default is white).
    """
    alignment = load_alignment(fasta_file)
    conservation_scores = calculate_conservation(alignment)

    plot_conserved_regions(conservation_scores, color, output_path_conserved)
