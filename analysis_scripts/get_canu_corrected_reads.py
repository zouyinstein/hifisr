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
import hifisr_functions.reads as hfreads
import hifisr_functions.reports as hfrps
import os



# Usage:
# Load the soft paths
soft_paths_file = sys.argv[1] 
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito
bait_ref = sys.argv[4] # mito_rotated_293434.fasta
sample_platform = sys.argv[5] # CLR, ONT
sample_num = int(sys.argv[6]) # 10000
threads = sys.argv[7] # 32


# create project folder by sample_index
if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists("reads"):
    os.makedirs("reads")
os.chdir("reads")


# correct reads by canu
hfreads.correct_reads_by_canu(sample_index, genome, bait_ref, sample_num, sample_platform, soft_paths_dict, threads, result_dir="canu_corrected_reads")
os.chdir("canu_corrected_reads")
id_length_qual_file, total_read_number, total_bases = hfrps.get_fastq_stats("corrected_" + genome, sample_index + "_" + genome + ".correctedReads.fasta.gz", soft_paths_dict, threads)
hfrps.plot_length_qual("corrected_" + genome, sample_platform, id_length_qual_file, total_read_number, total_bases)

os.chdir("../../..")