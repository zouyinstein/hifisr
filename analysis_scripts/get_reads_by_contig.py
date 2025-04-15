import sys
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
out_prefix = sys.argv[3] # mt_draft_contigs
contigs_absolute_path = sys.argv[4]
reads_absolute_path = sys.argv[5] # ATHiFi001.fastq
threads = sys.argv[6]

# check absolute paths
if not os.path.isabs(contigs_absolute_path):
    print("The path to the contigs reference fasta file is not an absolute path.")
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
if not os.path.exists(out_prefix):
    os.makedirs(out_prefix)
os.chdir(out_prefix)


# split the reads by contig
hfref.replace_fasta_id("contig", contigs_absolute_path, "contigs.fa")
command_1 = "ln -sf " + reads_absolute_path + " " + out_prefix + ".fastq"
ret = hfbase.get_cli_output_lines(command_1, side_effect = True)
count = hfreads.split_reads_by_contig(out_prefix + ".fastq", "HiFi", "contigs.fa", soft_paths_dict, threads)

# get the stats of the reads with plots
for i in range(1, int(count) + 1):
    id_length_qual_file, total_read_number, total_bases = hfrps.get_fastq_stats("contig_" + str(i), "contig_" + str(i) + ".fastq", soft_paths_dict, threads)
    hfrps.plot_length_qual("contig_" + str(i), "HiFi", id_length_qual_file, total_read_number, total_bases)

os.chdir("../../..")