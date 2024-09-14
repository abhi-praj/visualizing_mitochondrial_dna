import pandas as pd
from Bio import SeqIO
import plotly.express as px


def sequence_lengths(sequences: list[SeqIO.SeqRecord]) -> pd.DataFrame:
    """Extract sequence lengths from a list of sequences."""
    lengths = [len(seq.seq) for seq in sequences]
    ids = [seq.id for seq in sequences]
    return pd.DataFrame({'ID': ids, 'Length': lengths})


def plot_sequence_length_histogram(lengths_df: pd.DataFrame, output_path: str) -> None:
    """Plot and save a histogram of sequence lengths."""
    fig = px.histogram(lengths_df, x='Length', nbins=50, title="Sequence Length Distribution",
                       labels={'Length': 'Sequence Length'},
                       template='plotly_white')
    fig.update_layout(
        xaxis_title="Sequence Length",
        yaxis_title="Frequency",
        bargap=0.2
    )

    fig.write_html(output_path)


def sequence_length_distribution(fasta_file: str, output_image: str) -> None:
    """Load sequences from a FASTA file and generate a sequence length distribution plot."""
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    lengths_df = sequence_lengths(sequences)
    plot_sequence_length_histogram(lengths_df, output_image)
