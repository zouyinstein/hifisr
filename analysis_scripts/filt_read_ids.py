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
mito_id_length_qual_file = sys.argv[3] # mito_id_length_qual.txt
plastid_id_length_qual_file = sys.argv[4] # plastid_id_length_qual.txt
filt_length = int(sys.argv[5]) # 10000
filt_qual = float(sys.argv[6]) # 0

# create project folder by sample_index
if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists("reads"):
    os.makedirs("reads")
os.chdir("reads")

# filter the reads with length larger than 10K
id_length_qual_file_filt, filt_read_number, filt_bases = hfreads.filt_length_qual("L10K_mito", mito_id_length_qual_file, filt_length=filt_length, filt_qual=filt_qual)
hfrps.plot_length_qual("filt_L10K_mito", "HiFi", id_length_qual_file_filt, filt_read_number, filt_bases)
id_length_qual_file_filt, filt_read_number, filt_bases = hfreads.filt_length_qual("L10K_plastid", plastid_id_length_qual_file, filt_length=filt_length, filt_qual=filt_qual)
hfrps.plot_length_qual("filt_L10K_plastid", "HiFi", id_length_qual_file_filt, filt_read_number, filt_bases)

# move the filtered reads to filt_reads
if not os.path.exists("filt_reads"):
    os.makedirs("filt_reads")
command_1 = "mv filt*.* filt_reads"
hfbase.get_cli_output_lines(command_1, side_effect = True)

# extract the ids of the filtered reads
mito_ids_file = "filt_reads/filt_L10K_mito_ids.txt"
command_1 = "cut -f 1 filt_reads/filt_L10K_mito_id_length_qual.txt > " + mito_ids_file
hfbase.get_cli_output_lines(command_1, side_effect = True)
plastid_ids_file = "filt_reads/filt_L10K_plastid_ids.txt"
command_1 = "cut -f 1 filt_reads/filt_L10K_plastid_id_length_qual.txt > " + plastid_ids_file
hfbase.get_cli_output_lines(command_1, side_effect = True)

os.chdir("../..")