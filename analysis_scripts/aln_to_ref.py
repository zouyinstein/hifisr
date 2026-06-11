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
import os


# Usage: python
# Load the soft paths
soft_paths_file = sys.argv[1]
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito or plastid
node_prefix = sys.argv[4] # e.g., node_1
genome_absolute_path = sys.argv[5] # for bait genome
before_fasta_absolute_path = sys.argv[6] # draft.fasta



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
if not os.path.exists("aln_org"):
    os.makedirs("aln_org")
os.chdir("aln_org")
if not os.path.exists(genome):
    os.makedirs(genome)
os.chdir(genome)
if not os.path.exists(node_prefix):
    os.makedirs(node_prefix)
os.chdir(node_prefix)


# check the SVs by BLASTn
adj_count, rc_or_not, rot_or_not = hfref.aln_to_ref(genome, genome_absolute_path, before_fasta_absolute_path, genome + "_no_change.fasta", "all_sorted_blastn_alignments.txt", soft_paths_dict)
hfrps.convert_blastn_alignments_to_table("all_sorted_blastn_alignments.txt", genome + "_blastn_alignments.xlsx")
if adj_count == 0:
    print("No change.") # must be no change
    os.remove(genome + "_no_change.fasta")
    os.remove("all_sorted_blastn_alignments.txt")
else:
    print(f"Number of adjacencies detected: {adj_count}")
# if rc_or_not == 1:
#     print("The draft genome has been reverse complemented.")
# if rot_or_not == 1:
#     print("The draft genome has been rotated.")

os.chdir("../../../../..")