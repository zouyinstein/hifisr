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
import hifisr_functions.variants as hfvar
import hifisr_functions.references as hfref
import hifisr_functions.annotations as hfanno
import hifisr_functions.transfer as hftrans
import hifisr_functions.reports as hfrps
from Bio import SeqIO
import polars as pl
import pandas as pd
import subprocess
import sys
import os


# Usage: python
# Load the soft paths
soft_paths_file = sys.argv[1]
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
round_index = sys.argv[3] # run_1
chr_absolute_path = sys.argv[4]
mito_absolute_path = sys.argv[5]
plastid_absolute_path = sys.argv[6]


# check absolute paths
if not os.path.isabs(chr_absolute_path):
    print("The path to the nuclear reference fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(mito_absolute_path):
    print("The path to the mito fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(plastid_absolute_path):
    print("The path to the plastid fasta file is not an absolute path.")
    sys.exit(1)


# create project folder by sample_index
if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists("transfer_fragments"):
    os.makedirs("transfer_fragments")
os.chdir("transfer_fragments")
if not os.path.exists(round_index):
    os.makedirs(round_index)
os.chdir(round_index)


num_Chr, num_mito, num_plastid = hftrans.run_transfer_blastn(chr_absolute_path, mito_absolute_path, plastid_absolute_path, soft_paths_dict)
hftrans.merge_numt_nupt(num_Chr, num_mito, num_plastid)
hftrans.get_raw_transfer_groups(num_Chr)


os.chdir("../../..")
