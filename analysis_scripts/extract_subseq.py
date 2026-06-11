# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

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
