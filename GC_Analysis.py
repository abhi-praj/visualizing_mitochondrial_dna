import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.SeqUtils import gc_fraction
import numpy as np
from typing import Union


def display_gcc(sequences: Union[list, str]) -> None:
    """Display the GC content for various sequences."""
    if isinstance(sequences, list):
        gc_contents = [gc_fraction(seq.seq) for seq in sequences]
        gc_df = pd.DataFrame(
            {'ID': [seq.id for seq in sequences], 'GC_Content': gc_contents})
    else:
        sequences = [sequences]
        gc_contents = [gc_fraction(seq.seq) for seq in sequences]
        gc_df = pd.DataFrame(
            {'ID': [sequences[0].id], 'GC_Content': gc_contents})

    plt.figure(figsize=(10, 6))
    sns.histplot(gc_df['GC_Content'], bins=50, kde=True, color='skyblue',
                 edgecolor='black')
    plt.title('Distribution of GC Content in Mitochondrial DNA Sequences')
    plt.xlabel('GC Content Frequency')
    plt.ylabel('m-DNA Sequence Count')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()


def gc_skew(sequence: str, sequence_length: int) -> np.ndarray:
    """Calculate the GC skew for a given genome with a
     specified sequence length."""
    g_counts = np.array([sequence[i:i + sequence_length].count('G') for i in
                         range(len(sequence) - sequence_length + 1)])
    c_counts = np.array([sequence[i:i + sequence_length].count('C') for i in
                         range(len(sequence) - sequence_length + 1)])
    gc_skew_values = (g_counts - c_counts) / (g_counts + c_counts)
    return gc_skew_values


def plot_gc_skew(sequences: Union[list, str],
                 sequence_length: int = 1000) -> None:
    """Plot the GC skew for a list of genomes per each sequence input."""
    if isinstance(sequences, list):
        for seq_record in sequences:
            skew = gc_skew(str(seq_record.seq), sequence_length)
            plt.figure(figsize=(12, 6))
            plt.plot(skew, label=seq_record.id)
            plt.title(f'GC Skew for {seq_record.id}')
            plt.xlabel('Position in Sequence')
            plt.ylabel('GC Skew')
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.legend()
            plt.show()
    else:
        seq_record = sequences
        skew = gc_skew(str(seq_record.seq), sequence_length)
        plt.figure(figsize=(12, 6))
        plt.plot(skew, label=seq_record.id)
        plt.title(f'GC Skew for {seq_record.id}')
        plt.xlabel('Position in Sequence')
        plt.ylabel('GC Skew')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.show()


def plot_gc_skew_genome_set(sequences: list,
                            sequence_length: int = 1000) -> None:
    """Plot the GC skew for all genomes in a genome set as a
     single visualization set with a specific sequence length."""
    plt.figure(figsize=(14, 8))
    for seq_record in sequences:
        skew = gc_skew(str(seq_record.seq), sequence_length)
        plt.plot(skew, label=seq_record.id)

    plt.title('GC Skew for All Sequences')
    plt.xlabel('Position in Sequence')
    plt.ylabel('GC Skew')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1), ncol=1)
    plt.show()


def calculate_gc_content(sequence: str, sequence_length: int) -> np.ndarray:
    """Calculate the GC content for a given genome
     with a specified sequence len."""
    gc_content_values = np.array([
        gc_fraction(sequence[i:i + sequence_length]) for i in
        range(len(sequence) - sequence_length + 1)
    ])
    return gc_content_values


def gc_content_heatmap(sequences: Union[list, str], sequence_length: int = 1000) -> None:
    """Generate a heatmap of GC content for a list of genomes."""
    if isinstance(sequences, list):
        gc_content_matrix = []
        for seq_record in sequences:
            gc_content_values = calculate_gc_content(str(seq_record.seq), sequence_length)
            gc_content_matrix.append(gc_content_values)
        index_labels = [seq.id for seq in sequences]
    else:
        seq_record = sequences
        gc_content_matrix = [calculate_gc_content(str(seq_record.seq), sequence_length)]
        index_labels = [seq_record.id]

    gc_content_df = pd.DataFrame(gc_content_matrix, index=index_labels)

    plt.figure(figsize=(14, 8))
    sns.heatmap(gc_content_df, cmap='viridis', cbar=True)
    plt.title('GC Content Heatmap')
    plt.xlabel('Position in Sequence (Window Start)')
    plt.ylabel('Sequence ID')
    plt.show()


def cumulative_gc_content_plot(sequences: Union[list, str],
                               sequence_length: int = 1000) -> None:
    """Generate a cumulative plot of GC content for a list of genomes."""
    if isinstance(sequences, list):
        cumulative_gc = np.zeros(1)  # Initialize with a dummy array to be extended

        for seq_record in sequences:
            gc_content_values = calculate_gc_content(str(seq_record.seq), sequence_length)
            cumulative_gc = np.append(cumulative_gc, gc_content_values)

        cumulative_gc = cumulative_gc[1:]  # Remove the initial dummy array
    else:
        seq_record = sequences
        cumulative_gc = calculate_gc_content(str(seq_record.seq), sequence_length)

    plt.figure(figsize=(12, 6))
    plt.plot(np.arange(len(cumulative_gc)), cumulative_gc, color='skyblue', linewidth=2)
    plt.fill_between(np.arange(len(cumulative_gc)), cumulative_gc, color='skyblue', alpha=0.3)
    plt.title('Cumulative GC Content Plot')
    plt.xlabel('Position in Sequence')
    plt.ylabel('Cumulative GC Content')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.show()
