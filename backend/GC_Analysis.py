import pandas as pd
import plotly.express as px
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

    fig = px.histogram(gc_df, x='GC_Content', nbins=50, title="Distribution of GC Content in Mitochondrial DNA Genomes",
                       labels={'GC_Content': 'GC Content Frequency'},
                       template='plotly_white')
    fig.update_layout(
        xaxis_title="GC Content",
        yaxis_title="Genome Count",
        bargap=0.2
    )

    fig.write_html(output_image_path)
