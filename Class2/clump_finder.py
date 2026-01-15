import sys
from Bio import SeqIO
from collections import Counter

def find_clumps(sequence, k, L, t):
    """
    Detect (L, k, t)-clumps:
    k-mer appears at least t times within a window of length L.
    """
    clumps = set()
    n = len(sequence)

    for i in range(n - L + 1):
        window = sequence[i:i+L]
        counts = Counter()

        # Count k-mers in the window
        for j in range(L - k + 1):
            kmer = window[j:j+k]
            counts[kmer] += 1

        # Add k-mers that meet threshold t
        for kmer, count in counts.items():
            if count >= t:
                clumps.add(kmer)

    return clumps

def main():
    if len(sys.argv) < 2:
        print("Usage: python clump_finder.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]

    # Given parameters
    k = 8
    t = 3
    L = 1000

    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) == 0:
        print("No sequences found in FASTA file.")
        sys.exit(1)

    seq = str(records[0].seq).upper()

    print(f"Searching for ({L}, {k}, {t})-clumps...")

    clumps = find_clumps(seq, k, L, t)

    print(f"Number of clumps found: {len(clumps)}")
    print("Clumps:")
    for kmer in clumps:
        print(kmer)

if __name__ == "__main__":
    main()

