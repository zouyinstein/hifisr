import sys
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
genome_absolute_path = sys.argv[4] # for bait genome
before_fasta_absolute_path = sys.argv[5] # draft.fasta
reads_absolute_path = sys.argv[6] # ATHiFi001.fastq
threads = sys.argv[7]

# check absolute paths
if not os.path.isabs(genome_absolute_path):
    print("The path to the genome reference fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(before_fasta_absolute_path):
    print("The path to the draft fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(reads_absolute_path):
    print("The path to the reads fastq file is not an absolute path.")
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

hfref.replace_fasta_id(genome, before_fasta_absolute_path, "input_" + genome + ".fa")
rot_step = hfref.rotate_ref_to_non_repeat_region(genome, "input_" + genome + ".fa", soft_paths_dict, rotation=True)
hfbase.get_cli_output_lines("rm input_" + genome + ".fa")
hfref.flye_polish(genome, genome + "_rotated_" + str(rot_step) + ".fasta", genome + "_flye_polish", reads_absolute_path, soft_paths_dict, "HiFi", threads, correction=True)
# check the SVs by BLASTn
adj_count, rc_or_not, rot_or_not = hfref.aln_to_ref(genome, genome_absolute_path, genome + "_flye_polish.fasta", genome + "_flye_polish_aligned.fasta", "all_sorted_blastn_alignments.txt", soft_paths_dict)
hfrps.convert_blastn_alignments_to_table("all_sorted_blastn_alignments.txt", genome + "_flye_polish_aligned_blastn_alignments.xlsx")
if adj_count == 0:
    print("No change.")
if rc_or_not == 1:
    print("The draft genome has been reverse complemented.")
if rot_or_not == 1:
    print("The draft genome has been rotated.")

os.chdir("../../..")
