from hifisr_functions import hifisr_functions as hf
import concurrent.futures as cf
from Bio import SeqIO
import pandas as pd
import sys
import os


soft_paths_file = sys.argv[1]
soft_paths_dict = { line.split("\t")[0]:line.split("\t")[1] for line in hf.get_file_lines(soft_paths_file) }
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito, plastid
run_info = sys.argv[4] # run_1, run_2, etc.
genome_absolute_path = sys.argv[5]
reads_absolute_path = sys.argv[6]
threads = sys.argv[7]


if not os.path.exists(sample_index): os.mkdir(sample_index)
os.chdir(sample_index)
if not os.path.exists(genome): os.mkdir(genome)
os.chdir(genome)
if not os.path.exists(run_info): os.mkdir(run_info)
os.chdir(run_info)
current_dir = os.getcwd()

# run multi-threading blastn
if not os.path.exists("tmp_blastn_results"):
    os.mkdir("tmp_blastn_results")
command_1 = soft_paths_dict.get("seqkit") + " fq2fa " + reads_absolute_path + " > reads.fasta"
hf.get_cli_output_lines(command_1, side_effect = True)
hf.replace_reads_id("reads.fasta", "new_reads.fasta")
read_records = tuple(SeqIO.parse("new_reads.fasta", "fasta")) # fastq will be converted to fasta
read_ref_pairs = tuple([ (read_records[i], genome_absolute_path, soft_paths_dict) for i in range(len(read_records)) ]) 
if os.path.exists("all_sorted_blastn_alignments.txt"):
    os.remove("all_sorted_blastn_alignments.txt")
with cf.ThreadPoolExecutor(int(threads)) as tex:
    results = tex.map(hf.run_blastn_sorter, read_ref_pairs)
# merge all sorted blastn alignments
hf.get_cli_output_lines("rm -rf tmp_blastn_results", side_effect = True)
hf.get_type_and_subtype("all_sorted_blastn_alignments.txt", 5, out_dir="read_group_files") # split files by type and subtype
blastn_info_files = ["read_group_files/" + file for file in os.listdir("read_group_files") if file.startswith("type_")]
for blastn_info_file in blastn_info_files:
    hf.check_FL_and_multi(blastn_info_file, 5, out_dir="FL_read_group_files", id_dir="IDs", report_dir="Reports")
# summarize the results
command_1 = "wc -l FL_read_group_files/*_blastn_results.txt > all_FL_report.txt"
command_2 = "wc -l IDs/*_partial_ids.txt read_group_files/other_blastn_results.txt > all_rm_report.txt"
command_3 = "wc -l all_sorted_blastn_alignments.txt >> all_FL_report.txt"
commands = command_1 + ";" + command_2 + ";" + command_3
hf.get_cli_output_lines(commands, side_effect = True)
# get read groups in excel format
groups_1 = [["NA"]]
groups_2 = hf.get_next_groups(groups_1) # [["ins", "NA"], ["ref", "NA"], ["rep", "NA"]]
groups_3 = hf.get_next_groups(groups_2)
groups_4 = hf.get_next_groups(groups_3)
groups_5 = hf.get_next_groups(groups_4)
groups_list = [groups_1, groups_2, groups_3, groups_4, groups_5] # for easy manipulation
# analyze each subtype, and group the reads by (se1, ss2), (se1, ss2, se2, ss3), (se1, ss2, se2, ss3, se3, ss4), (se1, ss2, se2, ss3, se3, ss4, se4, ss5)
for num_align in range(2, 6): # 2, 3, 4, 5
    hf.get_subgroups(groups_list[num_align-1], num_align, input_dir="FL_read_group_files", out_dir="combined_excel", ID_subgroup_dir="ID_subgroup") # 1，2，3，4
if  not os.path.exists("backup_info"):
    os.mkdir("backup_info")
command_1 = "mv FL_read_group_files backup_info; mv IDs backup_info; mv Reports backup_info; mv read_group_files backup_info; mv all_sorted_blastn_alignments.txt backup_info; mv all_rm_report.txt backup_info; mv combined_excel/*_backup.xlsx backup_info"
command_2 = "find backup_info/FL_read_group_files -size 0 -delete"
command_3 = "find ID_subgroup -size 0 -delete"
commands = command_1 + ";" + command_2 + ";" + command_3
hf.get_cli_output_lines(commands, side_effect = True)
type_2_rep_file = current_dir + "/combined_excel/type_2_subtype_rep_NA_summary_anno.xlsx"
IDs_dir = current_dir + "/backup_info/IDs"
reads_fasta = current_dir + "/new_reads.fasta" # clean_mito.fasta, not fastq
ref_len = len(SeqIO.read(genome_absolute_path, "fasta").seq)
try: 
    df_ref = pd.read_excel(current_dir + "/combined_excel/type_2_subtype_ref_NA_summary_anno.xlsx", sheet_name='Sheet1')
    df_ref_filt = df_ref[(df_ref["(se1, ss2)"] == str((1, ref_len))) | (df_ref["(se1, ss2)"] == str((ref_len, 1)))]
    type_2_ref_index_list = [ str(i) for i in df_ref_filt["old_index"].to_list() ]
    type_2_ref_files = [ current_dir + "/ID_subgroup/type_2_subtype_ref_NA_subgroup_" + i + "_FL_ids.txt" for i in type_2_ref_index_list ]
except FileNotFoundError:
    type_2_ref_files = []
try:
    df_rep = pd.read_excel(current_dir + "/combined_excel/type_2_subtype_rep_NA_summary_anno.xlsx", sheet_name='Sheet1')
    df_rep_filt = df_rep[df_rep["mid_olp_1"] >= 1000]
    type_2_rep_large_index_list = [ str(i) for i in df_rep_filt["old_index"].to_list() ]
    type_2_rep_large_files = [ current_dir + "/ID_subgroup/type_2_subtype_rep_NA_subgroup_" + i + "_FL_ids.txt" for i in type_2_rep_large_index_list ]
except FileNotFoundError:
    type_2_rep_large_files = []
hf.plot_bubble_type_2_rep_raw(type_2_rep_file, IDs_dir, genome_absolute_path)
# use the first two for type_1_cov, use all three for variant calling
FL_ids_files = [IDs_dir + "/" + file for file in os.listdir(IDs_dir) if file.endswith("_FL_ids.txt")]
partial_ids_files = [IDs_dir + "/" + file for file in os.listdir(IDs_dir) if file.endswith("_partial_ids.txt")]
command_1 = "cat " + " ".join(FL_ids_files) + " > FL_ids.txt"
command_2 = "cat " + " ".join(partial_ids_files) + " > partial_ids.txt"
command_3 = soft_paths_dict.get("seqkit") + " grep -f FL_ids.txt " + reads_fasta + " > FL.fasta"
command_4 = soft_paths_dict.get("seqkit") + " grep -f partial_ids.txt " + reads_fasta + " > partial.fasta"
command_5 = "cut -f1 " + current_dir + "/backup_info/FL_read_group_files/type_1_subtype_NA_FL_blastn_results.txt " + " ".join(type_2_ref_files) + " " + " ".join(type_2_rep_large_files) + " > variant_cov_ids.txt"
command_6 = soft_paths_dict.get("seqkit") + " grep -f variant_cov_ids.txt " + reads_fasta + " > variant_cov.fasta"
commands = command_1 + " ; " + command_2 + " ; " + command_3 + " ; " + command_4 + " ; " + command_5 + " ; " + command_6
hf.get_cli_output_lines(commands, side_effect = True)
hf.cal_ID_coverage("FL", genome_absolute_path, "FL.fasta", soft_paths_dict, threads, tmp_dir="tmp_coverage")
hf.cal_ID_coverage("partial", genome_absolute_path, "partial.fasta", soft_paths_dict, threads, tmp_dir="tmp_coverage")
hf.cal_ID_coverage("variant", genome_absolute_path, "variant_cov.fasta", soft_paths_dict, threads, tmp_dir="tmp_coverage")
ref_len = len(SeqIO.read(genome_absolute_path, "fasta").seq)
hf.plot_coverage("FL_cov.txt", "partial_cov.txt", "variant_cov.txt", 1, ref_len)
df_short = hf.get_bcftools_frequency(genome_absolute_path, "variant_cov.fasta", soft_paths_dict, threads, tmp_dir="tmp_variants")
hf.get_pysam_frequency(df_short, genome_absolute_path, "variant_cov.fasta", soft_paths_dict, threads, top=100, tmp_dir="tmp_pysam")
os.chdir("../../..")
