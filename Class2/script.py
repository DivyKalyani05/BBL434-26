import sys
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt

# Sliding window generator
def sliding_windows(sequence, window_size, step):
    for start in range(0, len(sequence) - window_size + 1, step):
        end = start + window_size
        yield start, end, sequence[start:end]

# Count k-mers inside a sequence
def count_kmers(seq, k):
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        counts[kmer] += 1
    return counts

def main():
    # Ensure FASTA file is provided
    if len(sys.argv) < 2:
        print("Usage: python ori_kmer_enrichment.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]

    # Parameters for ORI signal detection
    k = 8
    window = 5000
    step = 500

    # Load genome sequence (first record only)
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) == 0:
        print("No sequences found in FASTA file.")
        sys.exit(1)

    sequence = str(records[0].seq).upper()

    positions = []
    max_kmer_counts = []

    # Scan genome using sliding windows
    for start, end, window_seq in sliding_windows(sequence, window, step):
        kmer_counts = count_kmers(window_seq, k)

        # Find the highest-occurring k-mer in this window
        max_count = max(kmer_counts.values())
        positions.append(start)
        max_kmer_counts.append(max_count)

    # Plot enrichment
    plt.figure(figsize=(12, 5))
    plt.plot(positions, max_kmer_counts)
    plt.xlabel("Genome Position (Start of Window)")
    plt.ylabel(f"Highest {k}-mer Count in Window")
    plt.title("ORI Signal Checker – K-mer Enrichment Profile")
    plt.grid(True)
    plt.tight_layout()

    # SAVE PLOT
    plt.savefig("class_2_plot.png", dpi=300)
    print("Plot saved as class_2_plot.png")

    # SHOW PLOT
    plt.show()

if __name__ == "__main__":
    main()
