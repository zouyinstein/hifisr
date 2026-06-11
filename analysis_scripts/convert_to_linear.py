# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

from Bio import SeqIO
import sys
import _bootstrap  # noqa: F401


prefix = sys.argv[1] # Amborella_linear_mito_v1
fasta_file_list = sys.argv[2:]

sequence_list = [ str(SeqIO.read(fasta_file, "fasta").seq) for fasta_file in fasta_file_list ]
with open(prefix + ".fasta", "wt") as fout:
    print(">", prefix, file=fout)
    print("".join(sequence_list), file=fout)
