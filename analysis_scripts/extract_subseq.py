import sys
import _bootstrap  # noqa: F401
import hifisr_functions.references as hfref
from Bio import SeqIO
import os


ref_fasta = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

ref_record = SeqIO.read(ref_fasta, "fasta")
subseq, _ = hfref.get_subseq(ref_record, start, end, flank=0)
print(">subseq")
print(subseq.seq)
