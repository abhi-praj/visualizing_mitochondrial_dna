from Bio import SeqIO
import GC_Analysis as GC_Analysis
import SequenceLengthDistribution as SeqLenDist
import CodonUsageAnalysis as CodonUsage
import data_validation as dv

# Replace sample_genome.fasta with your own fasta file
file_path = 'sample_genome.fasta'

# Validate file
valid = dv.validate_file(file_path)

if valid:
    def load_sequences(fasta_file: str):
        bit = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            bit.append(record)
        return bit
    sequences = load_sequences(file_path)

    # GC Content Analysis
    # GC_Analysis.display_gcc(sequences)
    # GC_Analysis.plot_gc_skew(sequences[0], sequence_length=1000, interval=2000)
    # GC_Analysis.plot_gc_skew_genome_set(sequences[:10], sequence_length=1000, interval=2000)
    # GC_Analysis.gc_content_heatmap(sequences[:10], sequence_length=1000, interval=1000)
    # GC_Analysis.cumulative_gc_content_plot(sequences[:10], sequence_length=1000)
    #
    # # Sequence Length Comparisons
    # lengths_df = SeqLenDist.sequence_lengths(sequences)
    # SeqLenDist.plot_sequence_length_histogram(lengths_df)
    # SeqLenDist.plot_sequence_length_boxplot(lengths_df)
    # SeqLenDist.plot_sequence_length_violin_plot(lengths_df)

    # Codon Usage Visualizations
    CodonUsage.analyze_codon_usage(sequences)

