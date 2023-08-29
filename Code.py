from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sequences = [
    SeqRecord(Seq("ACCTACGACGACGCATAGCAGCGAGCTAGCT"), id="Sequence1"),
    SeqRecord(Seq("AGCTACGACGACGCATAGCAGCGAGCTAGCT"), id="Sequence2"),
    SeqRecord(Seq("AGCTAGCTACGACGACTAATCGACGACGACATATATATAGAGTGACAGCT"), id="Sequence3"),
    # Add more sequences as needed
]
from Bio import pairwise2

# Calculate pairwise identity and store results in a dictionary
pairwise_identity = {}
for seq1 in sequences:
    for seq2 in sequences:
        alignment = pairwise2.align.globalxx(seq1.seq, seq2.seq, one_alignment_only=True)
        identity = alignment[0].score / len(seq1)
        pairwise_identity[(seq1.id, seq2.id)] = identity
import numpy as np

# Create an empty similarity matrix
similarity_matrix = np.zeros((len(sequences), len(sequences)))

# Fill the matrix with pairwise identities
for i, seq1 in enumerate(sequences):
    for j, seq2 in enumerate(sequences):
        similarity_matrix[i, j] = pairwise_identity[(seq1.id, seq2.id)]
import seaborn as sns
import matplotlib.pyplot as plt

# Create a heatmap
sns.set(font_scale=1)
sns.heatmap(similarity_matrix, annot=True, cmap="YlGnBu", xticklabels=[seq.id for seq in sequences], yticklabels=[seq.id for seq in sequences])

# Add labels and title
plt.xlabel("Sequences")
plt.ylabel("Sequences")
plt.title("Sequence Similarity Heatmap")

# Show the plot
plt.show()
