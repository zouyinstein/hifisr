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
import shutil


# Usage:
# Load the soft paths
soft_paths_file = sys.argv[1] 
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito
id_length_qual_file = sys.argv[4] # filt_L10K_mito_id_length_qual.txt, 
sample_num = int(sys.argv[5]) if len(sys.argv) > 5 else 10000
threads = sys.argv[6] if len(sys.argv) > 6 else "1"
sample_prefix = "sample_" + genome


# create project folder by sample_index
if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists("reads"):
    os.makedirs("reads")
os.chdir("reads")

# random sampling the reads with sample_num
id_length_qual_file_sampled, sample_read_number, sample_bases = hfreads.random_sampling(sample_prefix, id_length_qual_file, sample_number=sample_num)
hfrps.plot_length_qual(sample_prefix, "HiFi", id_length_qual_file_sampled, sample_read_number, sample_bases)

# move the sampled reads to sample_reads
if not os.path.exists("sample_reads"):
    os.makedirs("sample_reads")
for filename in [
    id_length_qual_file_sampled,
    sample_prefix + "_length_qual_distribution.pdf",
    sample_prefix + "_length_qual_2d_distribution.pdf",
]:
    if os.path.exists(filename):
        shutil.move(filename, os.path.join("sample_reads", filename))

# extract the ids of the sampled reads and extract the reads
ids_file = "sample_reads/" + sample_prefix + "_ids.txt"
sample_id_length_qual_file = "sample_reads/" + sample_prefix + "_id_length_qual.txt"
sample_fastq_file = "sample_reads/" + sample_prefix + ".fastq.gz"
sample_fastq_file_uncompressed = "sample_reads/" + sample_prefix + ".fastq"
if os.path.exists(sample_fastq_file_uncompressed):
    os.remove(sample_fastq_file_uncompressed)
fastq_candidates = [
    genome + ".fastq.gz",
    genome + ".fastq",
    sample_index + "_" + genome + ".fastq.gz",
    sample_index + "_" + genome + ".fastq",
]
fastq_file = next((path for path in fastq_candidates if os.path.exists(path)), None)
if fastq_file is None:
    print("Cannot find any genome reads fastq: " + ", ".join(fastq_candidates), file=sys.stderr)
    sys.exit(1)

pigz = soft_paths_dict.get("pigz", "pigz")
command_2 = "cut -f 1 " + sample_id_length_qual_file + " > " + ids_file
command_3 = soft_paths_dict.get("seqkit") + " grep -f " + ids_file + " " + fastq_file + " | " + pigz + " -p " + threads + " -c > " + sample_fastq_file
hfbase.get_cli_output_lines(command_2 + " && " + command_3, side_effect = True)

os.chdir("../..")
