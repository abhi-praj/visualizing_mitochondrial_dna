from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.SeqUtils import gc_fraction
import GCContentDistribution


def load_sequences(fasta_file):
    bit = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        bit.append(record)
    return bit


file_path = 'sample_genome.fasta'
sequences = load_sequences(file_path)


# Display the first few sequences
GCContentDistribution.display_gcc(sequences)
