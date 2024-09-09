import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO


def load_sequences(fasta_file: str):
    """Load sequences from a FASTA file."""
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(record)
    return sequences


def display_gcc(fasta_file: str, output_image_path: str) -> None:
    """Display the GC content for sequences in the given FASTA file."""
    sequences = load_sequences(fasta_file)
    gc_contents = [gc_fraction(seq.seq) for seq in sequences]
    gc_df = pd.DataFrame(
        {'ID': [seq.id for seq in sequences], 'GC_Content': gc_contents})

    plt.figure(figsize=(10, 6))
    sns.histplot(gc_df['GC_Content'], bins=50, kde=True, color='skyblue',
                 edgecolor='black')
    plt.title('Distribution of GC Content in Mitochondrial DNA Genomes')
    plt.xlabel('GC Content Frequency')
    plt.ylabel('m-DNA Genome Count')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.savefig(output_image_path)
    plt.close()
