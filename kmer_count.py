import sys
from Bio import SeqIO
from collections import Counter

# Check arguments
if len(sys.argv) < 3:
    print("Usage: python kmer_count.py <fasta_file> <k>")
    sys.exit(1)

fasta_file = sys.argv[1]
k = int(sys.argv[2])

kmer_counts = Counter()
total_kmers = 0

# Process each sequence in the FASTA file
for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq).upper()

    # Slide a window of length k
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmer_counts[kmer] += 1
        total_kmers += 1

# Print totals
print(f"Total k-mers counted (including repeats): {total_kmers}")
print(f"Total unique k-mers: {len(kmer_counts)}\n")

# Print each k-mer with counts
print(f"Unique {k}-mers and their counts:")
for kmer, count in kmer_counts.items():
    print(f"{kmer}: {count}")

