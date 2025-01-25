from hifisr_functions import hifisr_functions as hf
import sys
import os


soft_paths_file = sys.argv[1]
soft_paths_dict = { line.split("\t")[0]:line.split("\t")[1] for line in hf.get_file_lines(soft_paths_file) }
mito_absolute_path = sys.argv[2]
plastid_absolute_path = sys.argv[3]
reads_absolute_path = sys.argv[4] # ATHiFi001.fastq
sample_index = sys.argv[5] # ATHiFi001
threads = sys.argv[6]


if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists("reads"):
    os.makedirs("reads")
os.chdir("reads")
hf.replace_bait_id("mito", mito_absolute_path, "mito.fa")
hf.replace_bait_id("plastid", plastid_absolute_path, "plastid.fa")
command_1 = "ln -sf " + reads_absolute_path + " " + sample_index + ".fastq"
ret = hf.get_cli_output_lines(command_1, side_effect = True)
hf.split_mtpt_reads(sample_index, sample_index + ".fastq", "HiFi", "mito.fa", "plastid.fa", soft_paths_dict, threads)

length_list, qual_list, total_read_number, total_bases = hf.get_fastq_stats("all", sample_index + ".fastq", soft_paths_dict, threads)
hf.plot_length_qual("all", "HiFi", length_list, qual_list, total_read_number, total_bases)

filt_length = 1000
hf.filt_fastq_records(sample_index + "_mito_f1k", sample_index + "_mito.fastq", filt_length)
length_list, qual_list, total_read_number, total_bases = hf.get_fastq_stats("mito", sample_index + "_mito_f1k.fastq", soft_paths_dict, threads)
hf.plot_length_qual("mito", "HiFi", length_list, qual_list, total_read_number, total_bases)
hf.filt_fastq_records(sample_index + "_plastid_f1k", sample_index + "_plastid.fastq", filt_length)
length_list, qual_list, total_read_number, total_bases = hf.get_fastq_stats("plastid", sample_index + "_plastid_f1k.fastq", soft_paths_dict, threads)
hf.plot_length_qual("plastid", "HiFi", length_list, qual_list, total_read_number, total_bases)

command_1 = soft_paths_dict.get("seqkit") + " seq -ni " + sample_index + "_mito.fastq > " + sample_index + "_mito_fastq_ids.txt"
command_2 = "rm -rf " + sample_index + "_mito.fastq"
command_3 = soft_paths_dict.get("seqkit") + " seq -ni " + sample_index + "_plastid.fastq > " + sample_index + "_plastid_fastq_ids.txt"
command_4 = "rm -rf " + sample_index + "_plastid.fastq"
commands = command_1 + " ; " + command_2 + " ; " + command_3 + " ; " + command_4
ret = hf.get_cli_output_lines(commands, side_effect = True)
os.chdir("../..")