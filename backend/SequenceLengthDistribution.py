import pandas as pd
from Bio import SeqIO
import seaborn as sns
from matplotlib import pyplot as plt


def sequence_lengths(sequences: list[SeqIO.SeqRecord]) -> pd.DataFrame:
    """
    Extract sequence lengths from a list of sequences.
    """
    lengths = [len(seq.seq) for seq in sequences]
    ids = [seq.id for seq in sequences]
    return pd.DataFrame({'ID': ids, 'Length': lengths})


def plot_sequence_length_histogram(lengths_df: pd.DataFrame,
                                   output_path: str) -> None:
    """
    Plot and save a histogram of sequence lengths.
    """
    plt.figure(figsize=(10, 6))
    sns.histplot(lengths_df['Length'], bins=50, kde=True, color='skyblue',
                 edgecolor='black')
    plt.title('Sequence Length Distribution')
    plt.xlabel('Sequence Length')
    plt.ylabel('Frequency')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig(output_path)
    plt.close()


def sequence_length_distribution(fasta_file: str, output_image: str) -> None:
    """
    Load sequences from a fasta file and generate a sequence length distribution plot.
    """
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    lengths_df = sequence_lengths(sequences)
    plot_sequence_length_histogram(lengths_df, output_image)
