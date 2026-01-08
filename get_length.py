import sys
from Bio import SeqIO

# Check if the user gave a file name
if len(sys.argv) < 2:
    print("Usage: python get_length.py <fasta_file>")
    sys.exit(1)

fasta_file = sys.argv[1]

# Read and process sequences
for record in SeqIO.parse(fasta_file, "fasta"):
    print(f"ID: {record.id}")
    print(f"Length: {len(record.seq)}")
    print()

