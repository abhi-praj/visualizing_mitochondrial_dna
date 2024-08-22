import pandas as pd
from Bio.SeqUtils import CodonAdaptationIndex
import matplotlib.pyplot as plt
import seaborn as sns


def extract_codon_usage(sequence) -> pd.Series:
    """Extract codon usage from a sequence and return it as a Pandas Series."""
    sequence = ''.join(c for c in sequence.seq if c in 'ATCG')

    length = len(sequence)

    if length % 3 != 0:
        new_length = length - (length % 3)
        sequence = sequence[:new_length]
    codon_usage = CodonAdaptationIndex([sequence])

    return pd.Series(codon_usage)


def create_codon_usage_df(sequences) -> pd.DataFrame:
    """Create a DataFrame of codon usage from a list of sequences."""
    codon_usage_data = []
    for sequence in sequences:
        codon_usage_series = extract_codon_usage(sequence)
        codon_usage_data.append(codon_usage_series)
    return pd.DataFrame(codon_usage_data)


def plot_codon_usage_histogram(codon_usage_df: pd.DataFrame) -> None:
    """Plot a histogram of codon usage."""
    plt.figure(figsize=(12, 6))
    codon_usage_df.plot(kind='hist', bins=50, edgecolor='black', alpha=0.7)
    plt.title('Histogram of Codon Usage')
    plt.xlabel('Codon Usage')
    plt.ylabel('Frequency')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()


def plot_codon_usage_boxplot(codon_usage_df: pd.DataFrame) -> None:
    """Plot a boxplot of codon usage."""
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=codon_usage_df, color='skyblue')
    plt.title('Box Plot of Codon Usage')
    plt.xlabel('Codon Usage')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()


def plot_codon_usage_violin_plot(codon_usage_df: pd.DataFrame) -> None:
    """Plot a violin plot of codon usage."""
    plt.figure(figsize=(12, 6))
    sns.violinplot(data=codon_usage_df, color='skyblue')
    plt.title('Violin Plot of Codon Usage')
    plt.xlabel('Codon Usage')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()


def analyze_codon_usage(sequences: list) -> None:
    """Analyze and plot codon usage from a FASTA file."""
    codon_usage_df = create_codon_usage_df(sequences)

    # Plot the codon usage
    plot_codon_usage_histogram(codon_usage_df)
    plot_codon_usage_boxplot(codon_usage_df)
    plot_codon_usage_violin_plot(codon_usage_df)
