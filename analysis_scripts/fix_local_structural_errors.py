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
import hifisr_functions.reads as hfreads
import hifisr_functions.variants as hfvar
import hifisr_functions.annotations as hfanno
import hifisr_functions.reports as hfrps
from Bio import SeqIO
import os


# Usage: python get_variants_in_reads.py soft_paths.txt ATHiFi001 mito run_1 /path/to/ref/mito_rotated.fasta /path/to/reads/sample_4000_mito.fastq 32
# Load the soft paths
soft_paths_file = sys.argv[1] 
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito
run_info = sys.argv[4] # run_1
genome_absolute_path = sys.argv[5]
reads_absolute_path = sys.argv[6]
threads = sys.argv[7]
local_start = sys.argv[8]
local_end = sys.argv[9]


# check absolute paths
if not os.path.isabs(genome_absolute_path):
    print("The path to the genome reference fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(reads_absolute_path):
    print("The path to the reads fastq file is not an absolute path.")
    sys.exit(1)

# create project folder by sample_index
if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists(genome):
    os.makedirs(genome)
os.chdir(genome)
if not os.path.exists(run_info):
    os.makedirs(run_info)
os.chdir(run_info)
current_dir = os.getcwd()


hfrps.plot_coverage("FL_cov.txt", "partial_cov.txt", "variant_cov.txt", int(local_start), int(local_end), fig_length=12, fig_height=3)


os.chdir("../../..")
