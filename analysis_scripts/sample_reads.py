import sys
sys.path.append('/mnt/software/bio/hifisr/dev') # Add the path to the hifisr_functions.py file
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
id_length_qual_file = sys.argv[4] # filt_L10K_mito_id_length_qual.txt, 
sample_num = int(sys.argv[5]) # 4000


# create project folder by sample_index
if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists("reads"):
    os.makedirs("reads")
os.chdir("reads")

# random sampling the reads with sample_num
id_length_qual_file_sampled, sample_read_number, sample_bases = hfreads.random_sampling("sample_" + str(sample_num) + "_" + genome, id_length_qual_file, sample_number=sample_num)
hfrps.plot_length_qual("sample_" + str(sample_num) + "_" + genome, "HiFi", id_length_qual_file_sampled, sample_read_number, sample_bases)

# move the sampled reads to sample_reads
if not os.path.exists("sample_reads"):
    os.makedirs("sample_reads")
command_1 = "mv sample*.* sample_reads"
hfbase.get_cli_output_lines(command_1, side_effect = True)

# extract the ids of the sampled reads and extract the reads
ids_file = "sample_reads/sample_" + str(sample_num) + "_" + genome + "_ids.txt"
fastq_file = sample_index + "_" + genome + ".fastq"
command_1 = "cut -f 1 sample_reads/sample_" + str(sample_num) + "_" + genome + "_id_length_qual.txt > " + ids_file
command_2 = soft_paths_dict.get("seqkit") + " grep -f " + ids_file + " " + fastq_file + " > sample_reads/sample_" + str(sample_num) + "_" + genome + ".fastq"
hfbase.get_cli_output_lines(command_1 + " && " + command_2, side_effect = True)

os.chdir("../..")