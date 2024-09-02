from Bio import AlignIO, SeqIO
from Bio.Phylo import draw
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord


def preprocess_sequences(fasta_file: str, output_file: str) -> None:
    """
    Preprocess sequences to ensure they are of the same length.
    """
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
    """
    Load a multiple sequence alignment from a FASTA file.
    """
    temp_fasta_file = "preprocessed_phylogenetic" + fasta_file
    preprocess_sequences(fasta_file, temp_fasta_file)

    alignment = AlignIO.read(temp_fasta_file, "fasta")
    return alignment


def construct_phylogenetic_tree(alignment: MultipleSeqAlignment) -> Phylo:
    """
    Construct a phylogenetic tree from an alignment.
    """
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)

    return tree


def plot_phylogenetic_tree(tree: Phylo) -> None:
    """
    Plot the phylogenetic tree.
    """
    plt.figure(figsize=(10, 8))
    draw(tree, do_show=False)
    plt.title("Phylogenetic Tree")
    plt.show()


def phylogenetic_tree_pipeline(fasta_file: str) -> None:
    """
    Load sequences, construct a phylogenetic tree, and plot it.
    """
    alignment = load_alignment(fasta_file)
    tree = construct_phylogenetic_tree(alignment)
    plot_phylogenetic_tree(tree)
