import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

codon_to_amino_acid = {
    'GCT': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
    'CGT': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
    'AGA': 'Arginine', 'AGG': 'Arginine',
    'AAT': 'Asparagine', 'AAC': 'Asparagine',
    'GAT': 'Aspartic Acid', 'GAC': 'Aspartic Acid',
    'TGT': 'Cysteine', 'TGC': 'Cysteine',
    'CAA': 'Glutamine', 'CAG': 'Glutamine',
    'GAA': 'Glutamic Acid', 'GAG': 'Glutamic Acid',
    'GGT': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine',
    'CAT': 'Histidine', 'CAC': 'Histidine',
    'ATT': 'Isoleucine', 'ATC': 'Isoleucine', 'ATA': 'Isoleucine',
    'CTT': 'Leucine', 'CTC': 'Leucine', 'CTA': 'Leucine', 'CTG': 'Leucine',
    'TTA': 'Leucine', 'TTG': 'Leucine',
    'AAA': 'Lysine', 'AAG': 'Lysine',
    'ATG': 'Methionine',
    'TTT': 'Phenylalanine', 'TTC': 'Phenylalanine',
    'CCT': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
    'TCT': 'Serine', 'TCC': 'Serine', 'TCA': 'Serine', 'TCG': 'Serine',
    'AGT': 'Serine', 'AGC': 'Serine',
    'ACT': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine',
    'ACG': 'Threonine',
    'TGG': 'Tryptophan',
    'TAT': 'Tyrosine', 'TAC': 'Tyrosine',
    'GTT': 'Valine', 'GTC': 'Valine', 'GTA': 'Valine', 'GTG': 'Valine',
    'TAA': 'STOP', 'TAG': 'STOP', 'TGA': 'STOP'
}


def extract_codon_usage(sequence) -> pd.Series:
    """
    Extract codon usage from a sequence and return it as a Pandas Series.
    """
    sequence = ''.join(c for c in sequence.seq if c in 'ATCG')

    length = len(sequence)

    if length % 3 != 0:
        new_length = length - (length % 3)
        sequence = sequence[:new_length]

    codon_counts = {}
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        if codon in codon_to_amino_acid:
            amino_acid = codon_to_amino_acid[codon]
            if amino_acid in codon_counts:
                codon_counts[amino_acid] += 1
            else:
                codon_counts[codon_to_amino_acid[codon]] = 1

    return pd.Series(codon_counts)


def create_codon_usage_df(sequences) -> pd.DataFrame:
    """
    Create a DataFrame of codon usage from a list of sequences.
    """
    codon_usage_data = []
    for sequence in sequences:
        codon_usage_series = extract_codon_usage(sequence)
        codon_usage_data.append(codon_usage_series)
    return pd.DataFrame(codon_usage_data)


def plot_codon_usage_histogram(codon_usage_df: pd.DataFrame, output_path: str) -> None:
    """
    Plot a histogram of codon usage and save it as an image file.
    """
    codon_usage_df.fillna(0, inplace=True)
    plt.figure(figsize=(12, 6))
    codon_usage_df.plot(kind='hist', bins=50, edgecolor='black', alpha=0.7)
    plt.title('Histogram of Amino Acid Frequencies in Sequences')
    plt.xlabel('Number of Occurrences of Amino Acids per Sequence')
    plt.ylabel('Number of Sequences')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Save the plot as an image file
    plt.savefig(output_path)
    plt.close()  # Close the plot to free up memory


def analyze_codon_usage(filepath: str, output_path: str) -> None:
    """
    Analyze codon usage from a FASTA file and save the plot as an image.
    """
    sequences = list(SeqIO.parse(filepath, "fasta"))
    codon_usage_df = create_codon_usage_df(sequences)

    # Save the codon usage histogram plot as an image file
    plot_codon_usage_histogram(codon_usage_df, output_path)
