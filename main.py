from Bio import SeqIO
import GC_Analysis



def load_sequences(fasta_file):
    bit = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        bit.append(record)
    return bit


# Replace sample_genome.fasta with your own fasta file
file_path = 'sample_genome.fasta'
sequences = load_sequences(file_path)

# Display the GC Content Distribution the fasta file
# GC_Analysis.display_gcc(sequences)

# sequence_length represents the length of each sequence of a genome that is analyzed for GC
# Lower the sequence length, the more accurate but nosier the visualization is
# Default sequence length is 1000

# Plot the GC skew across various genomes in separate visualizations
# GC_Analysis.plot_gc_skew(sequences[0], sequence_length=1000, interval=2000)

# Plot the GC skew across various genomes in one big visualization
# GC_Analysis.plot_gc_skew_genome_set(sequences[:10], sequence_length=1000, interval=2000)

# Plot a heatmap of the GC content across various sequences of the genome for various genomes.
GC_Analysis.gc_content_heatmap(sequences[:10], sequence_length=1000, interval=1000)

# Plot Cumulative GC Content across all sequences
# GC_Analysis.cumulative_gc_content_plot(sequences[:10], sequence_length=1000)
