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
import hifisr_functions.base as hfbase
import hifisr_functions.references as hfref
import hifisr_functions.reports as hfrps
import sys
import os


# Usage: python
# Load the soft paths
soft_paths_file = sys.argv[1]
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito or plastid
genome_absolute_path = sys.argv[4] # for bait genome
before_fasta_absolute_path = sys.argv[5] # draft.fasta
flipped_fasta_path = sys.argv[6]
flip_start = sys.argv[7]
flip_end = sys.argv[8]

# check absolute paths
if not os.path.isabs(genome_absolute_path):
    print("The path to the genome reference fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(before_fasta_absolute_path):
    print("The path to the draft fasta file is not an absolute path.")
    sys.exit(1)

# create project folder by sample_index
if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists("draft_assembly"):
    os.makedirs("draft_assembly")
os.chdir("draft_assembly")
if not os.path.exists(genome):
    os.makedirs(genome)
os.chdir(genome)

# flipping via perfect repeats
hfref.get_flipped_fasta(before_fasta_absolute_path, flipped_fasta_path, flip_start, flip_end)
adj_count, rc_or_not, rot_or_not = hfref.aln_to_ref(genome, genome_absolute_path, flipped_fasta_path, genome + "_flye_polish_aligned_flipped_IR_aligned.fasta", "all_sorted_blastn_alignments.txt", soft_paths_dict)
hfrps.convert_blastn_alignments_to_table("all_sorted_blastn_alignments.txt", genome + "_flye_polish_aligned_flipped_IR_aligned_blastn_alignments.xlsx") # Overwrite the previous file
if adj_count == 0:
    print("No change.")
if rc_or_not == 1:
    print("The draft genome has been reverse complemented.")
if rot_or_not == 1:
    print("The draft genome has been rotated.")


os.chdir("../../..")
