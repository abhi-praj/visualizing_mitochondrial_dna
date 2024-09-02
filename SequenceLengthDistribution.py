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


def plot_sequence_length_histogram(lengths_df: pd.DataFrame) -> None:
    """
    Plot a histogram of sequence lengths.
    """
    plt.figure(figsize=(10, 6))
    sns.histplot(lengths_df['Length'], bins=50, kde=True, color='skyblue', edgecolor='black')
    plt.title('Sequence Length Distribution')
    plt.xlabel('Sequence Length')
    plt.ylabel('Frequency')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()


def plot_sequence_length_boxplot(lengths_df: pd.DataFrame) -> None:
    """
    Plot a boxplot of sequence lengths.
    """
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=lengths_df['Length'], color='skyblue')
    plt.title('Box Plot of Sequence Lengths')
    plt.xlabel('Sequence Length')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()


def plot_sequence_length_violin_plot(lengths_df: pd.DataFrame) -> None:
    """
    Plot a violin plot of sequence lengths.
    """
    plt.figure(figsize=(10, 6))
    sns.violinplot(x=lengths_df['Length'], color='skyblue')
    plt.title('Violin Plot of Sequence Lengths')
    plt.xlabel('Sequence Length')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

