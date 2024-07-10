from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.SeqUtils import gc_fraction


def display_gcc(sequence: list) -> None:
    for seq in sequence[:30]:
        print(f"ID: {seq.id}")
        print(f"Sequence: {seq.seq[:50]}...")
        print()

    gc_contents = [gc_fraction(seq.seq) for seq in sequence]

    gc_df = pd.DataFrame(
        {'ID': [seq.id for seq in sequence], 'GC_Content': gc_contents})

    print(gc_df.head())

    plt.figure(figsize=(10, 6))
    sns.histplot(gc_df['GC_Content'], bins=50, kde=True, color='skyblue',
                 edgecolor='black')
    plt.title('Distribution of GC Content in Mitochondrial DNA Sequences')
    plt.xlabel('GC Content Frequency')
    plt.ylabel('Mitochondrial DNA Sequence Count')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()
