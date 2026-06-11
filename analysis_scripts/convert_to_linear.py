from Bio import SeqIO
import sys
import _bootstrap  # noqa: F401


prefix = sys.argv[1] # Amborella_linear_mito_v1
fasta_file_list = sys.argv[2:]

sequence_list = [ str(SeqIO.read(fasta_file, "fasta").seq) for fasta_file in fasta_file_list ]
with open(prefix + ".fasta", "wt") as fout:
    print(">", prefix, file=fout)
    print("".join(sequence_list), file=fout)
