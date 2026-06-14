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
import hifisr_functions.reports as hfrps
import os


# Usage: python get_mtpt_reads.py soft_paths_file sample_index mito_absolute_path plastid_absolute_path reads_absolute_path threads
# Load the soft paths
soft_paths_file = sys.argv[1]
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
mito_absolute_path = sys.argv[3]
plastid_absolute_path = sys.argv[4]
reads_absolute_path = sys.argv[5] # ATHiFi001.fastq
threads = sys.argv[6]

# check absolute paths
if not os.path.isabs(mito_absolute_path):
    print("The path to the mito reference fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(plastid_absolute_path):
    print("The path to the plastid reference fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(reads_absolute_path):
    print("The path to the reads fastq file is not an absolute path.")
    sys.exit(1)

# create project folder by sample_index
if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists("reads"):
    os.makedirs("reads")
os.chdir("reads")

# split the reads into mito and plastid reads
hfref.replace_fasta_id("mito", mito_absolute_path, "mito.fa")
hfref.replace_fasta_id("plastid", plastid_absolute_path, "plastid.fa")
split_prefix = "split_mtpt"
if reads_absolute_path.endswith(".gz"):
    raw_reads_link = "all.fastq.gz"
else:
    raw_reads_link = "all.fastq"
cleanup_files = [
    "all.fastq", "all.fastq.gz",
    "mito.fastq", "mito.fastq.gz",
    "plastid.fastq", "plastid.fastq.gz",
    split_prefix + "_mito.fastq", split_prefix + "_plastid.fastq",
    sample_index + ".fastq", sample_index + ".fastq.gz",
    sample_index + "_mito.fastq", sample_index + "_mito.fastq.gz",
    sample_index + "_plastid.fastq", sample_index + "_plastid.fastq.gz",
]
hfbase.run_checked("rm -f " + " ".join(cleanup_files))
command_1 = "ln -sf " + reads_absolute_path + " " + raw_reads_link
hfbase.get_cli_output_lines(command_1, side_effect = True)
hfreads.split_mtpt_reads(split_prefix, raw_reads_link, "HiFi", "mito.fa", "plastid.fa", soft_paths_dict, threads)
os.replace(split_prefix + "_mito.fastq", "mito.fastq")
os.replace(split_prefix + "_plastid.fastq", "plastid.fastq")

id_length_qual_file, total_read_number, total_bases = hfrps.get_fastq_stats("all", raw_reads_link, soft_paths_dict, threads)
hfrps.plot_length_qual("all", "HiFi", id_length_qual_file, total_read_number, total_bases)
id_length_qual_file, total_read_number, total_bases = hfrps.get_fastq_stats("mito", "mito.fastq", soft_paths_dict, threads)
hfrps.plot_length_qual("mito", "HiFi", id_length_qual_file, total_read_number, total_bases)
id_length_qual_file, total_read_number, total_bases = hfrps.get_fastq_stats("plastid", "plastid.fastq", soft_paths_dict, threads)
hfrps.plot_length_qual("plastid", "HiFi", id_length_qual_file, total_read_number, total_bases)

pigz = soft_paths_dict.get("pigz", "pigz")
command_1 = pigz + " -p " + threads + " -f mito.fastq plastid.fastq"
hfbase.get_cli_output_lines(command_1, side_effect = True)

os.chdir("../..")
