import sys
sys.path.append('/mnt/software/bio/hifisr/dev') # Add the path to the hifisr_functions.py file
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
if reads_absolute_path.endswith(".gz"):
    command_1 = "ln -sf " + reads_absolute_path + " " + sample_index + ".fastq.gz"
    ret = hfbase.get_cli_output_lines(command_1, side_effect = True)
    hfreads.split_mtpt_reads(sample_index, sample_index + ".fastq.gz", "HiFi", "mito.fa", "plastid.fa", soft_paths_dict, threads)
    # get the stats of the reads with plots
    id_length_qual_file, total_read_number, total_bases = hfrps.get_fastq_stats("all", sample_index + ".fastq.gz", soft_paths_dict, threads)
    hfrps.plot_length_qual("all", "HiFi", id_length_qual_file, total_read_number, total_bases)
else:
    command_1 = "ln -sf " + reads_absolute_path + " " + sample_index + ".fastq"
    ret = hfbase.get_cli_output_lines(command_1, side_effect = True)
    hfreads.split_mtpt_reads(sample_index, sample_index + ".fastq", "HiFi", "mito.fa", "plastid.fa", soft_paths_dict, threads)
    # get the stats of the reads with plots
    id_length_qual_file, total_read_number, total_bases = hfrps.get_fastq_stats("all", sample_index + ".fastq", soft_paths_dict, threads)
    hfrps.plot_length_qual("all", "HiFi", id_length_qual_file, total_read_number, total_bases)
    
id_length_qual_file, total_read_number, total_bases = hfrps.get_fastq_stats("mito", sample_index + "_mito.fastq", soft_paths_dict, threads)
hfrps.plot_length_qual("mito", "HiFi", id_length_qual_file, total_read_number, total_bases)
id_length_qual_file, total_read_number, total_bases = hfrps.get_fastq_stats("plastid", sample_index + "_plastid.fastq", soft_paths_dict, threads)
hfrps.plot_length_qual("plastid", "HiFi", id_length_qual_file, total_read_number, total_bases)

os.chdir("../..")