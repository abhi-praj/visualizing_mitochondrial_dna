import os
import plotly.graph_objects as go
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
from Bio import SeqIO
from typing import Tuple


def preprocess_sequences(fasta_file: str, output_file: str) -> None:
    """Preprocess sequences to ensure they are of the same length."""
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
    """Load a multiple sequence alignment from a FASTA file."""
    temp_fasta_file = "preprocessed_conserved_" + os.path.basename(fasta_file)
    preprocess_sequences(fasta_file, temp_fasta_file)

    alignment = AlignIO.read(temp_fasta_file, "fasta")
    os.remove(temp_fasta_file)
    return alignment


def plot_conserved_regions(alignment: MultipleSeqAlignment, output_image_path: str) -> None:
    """Plot conserved regions in an alignment and save as an interactive Plotly plot."""
    alignment_array = np.array([list(record) for record in alignment], np.character)
    num_sequences, alignment_length = alignment_array.shape

    conservation = []
    for i in range(alignment_length):
        column = alignment_array[:, i]
        unique, counts = np.unique(column, return_counts=True)
        most_common_base = unique[np.argmax(counts)]
        conservation.append(np.sum(column == most_common_base) / num_sequences)

    fig = go.Figure(data=go.Scatter(x=np.arange(1, alignment_length + 1), y=conservation, mode='lines',
                                    line=dict(color='royalblue', width=2)))
    fig.update_layout(
        title="Conservation of Genomic Regions",
        xaxis_title="Position in Alignment",
        yaxis_title="Proportion of Conservation",
        template='plotly_white'
    )

    fig.write_html(output_image_path)


def conserved_regions(fasta_file: str, output_image_path: str) -> None:
    """Display conserved regions for sequences in a FASTA file."""
    alignment = load_alignment(fasta_file)
    plot_conserved_regions(alignment, output_image_path)
