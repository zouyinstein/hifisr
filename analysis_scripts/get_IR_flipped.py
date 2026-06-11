# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
import sys
import _bootstrap  # noqa: F401


input_fasta = sys.argv[1]
if len(sys.argv) != 4 and len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <input.fasta> <start> <end>", file=sys.stderr)
    sys.exit(1)

if len(sys.argv) == 2:
    command = "blastn -query " + input_fasta + " -subject " + input_fasta + " -outfmt 6 | head "
    subprocess.run(command, shell=True)
    exit(0)

if len(sys.argv) == 4:
    flip_start = int(sys.argv[2])
    flip_end = int(sys.argv[3])
    input_record = SeqIO.read(input_fasta, "fasta")
    input_id = input_record.id
    output_id = input_record.id + "_flipped_IR_" + str(flip_start) + "_" + str(flip_end)
    output_seq = (
        input_record.seq[0:flip_start-1]
        + input_record.seq[flip_start-1:flip_end].reverse_complement()
        + input_record.seq[flip_end:len(input_record.seq)]
    )
    output_record = SeqRecord(output_seq, id=output_id, description="")
    SeqIO.write(output_record, sys.stdout, "fasta")
