# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

import hifisr_functions.base as hfbase
import hifisr_functions.reads as hfreads
from Bio import SeqIO
import pandas as pd
import polars as pl
from collections import OrderedDict
import concurrent.futures as cf
import math
import sys
import os
import tempfile

# Function purity marker. "pure" means deterministic from explicit inputs with
# no file, shell, environment, logging, or input-mutation side effects.
FUNCTION_PURITY = {
    "get_tmp_root": "impure",
    "Index_label_alignments": "impure",
    "run_blastn_sorter_single": "impure",
    "run_multi_threads_blastn": "impure",
    "get_type_and_subtype": "impure",
    "check_FL_and_multi": "impure",
    "get_next_groups": "pure",
    "match_se1_ss2": "impure",
    "match_se1_ss2_se2_ss3": "impure",
    "match_se1_ss2_se2_ss3_se3_ss4": "impure",
    "match_se1_ss2_se2_ss3_se3_ss4_se4_ss5": "impure",
    "get_subgroups": "impure",
    "summarize_blastn_results": "impure",
    "get_cov_reads": "impure",
    "run_bcftools": "impure",
    "run_multi_threads_bcftools": "impure",
    "snv_or_indel": "impure",
}


def get_tmp_root():
    return os.environ.get("HIFISR_TMPDIR", tempfile.gettempdir())


class Index_label_alignments():
    def __init__(self, blastn_alignments_lines):
        self.blastn_alignments_lines = blastn_alignments_lines
        self.blastn_index_label_alignments = self.index_label_alignments()
    
    def index_label_alignments(self):
        if len(self.blastn_alignments_lines) == 0:
            return []
        blastn_index_label_alignments = [ [i, 1, ";".join(self.blastn_alignments_lines[i].split("\t"))] for i in range(len(self.blastn_alignments_lines))]
        return blastn_index_label_alignments

    def get_sorted_one_alignments(self):
        one_alignments = []
        for i in range(len(self.blastn_index_label_alignments)):
            if self.blastn_index_label_alignments[i][1] == 1:
                one_alignments.append(self.blastn_index_label_alignments[i])
        one_alignments.sort(key=lambda x: int(x[2].split(";")[6]))
        return one_alignments

    def pop_contained_alignments(self):
        count = 0
        for i in range(len(self.blastn_index_label_alignments)):
            if self.blastn_index_label_alignments[i][1] == 1:
                i_start = int(self.blastn_index_label_alignments[i][2].split(";")[6])
                i_end = int(self.blastn_index_label_alignments[i][2].split(";")[7])
                for j in range(len(self.blastn_index_label_alignments)):
                    if j == i:
                        continue
                    if self.blastn_index_label_alignments[j][1] == 1:
                        j_start = int(self.blastn_index_label_alignments[j][2].split(";")[6])
                        j_end = int(self.blastn_index_label_alignments[j][2].split(";")[7])
                        if i_start <= j_start and j_end <= i_end:
                            self.blastn_index_label_alignments[j][1] = 0
                            count += 1
                        elif j_start <= i_start and i_end <= j_end:
                            self.blastn_index_label_alignments[i][1] = 0
                            count += 1
        return count


def run_blastn_sorter_single(read_record, ref_fasta, results_filename, soft_paths_dict, tmp_blastn_dir):
    with open(tmp_blastn_dir + "/" + read_record.id + ".fasta", "wt") as fout:
        print(read_record.format("fasta"), file=fout)
    q_len = len(read_record.seq)
    ref_records = list(SeqIO.parse(ref_fasta, "fasta")) # ref_length = len(ref_records[0].seq)
    if len(ref_records) > 1:
        print("Warning: multiple records in reference file", file=sys.stderr)
    command_1 = soft_paths_dict.get("blastn") + " -query " + tmp_blastn_dir + "/" + read_record.id + ".fasta -subject " + ref_fasta + " -outfmt 6"
    blastn_alignments_lines = hfbase.get_cli_output_lines(command_1, side_effect = False)
    if len(blastn_alignments_lines) == 0:
        with open("reads_with_no_alignments.txt", "at") as fout:
            print(read_record.id, file=fout)
    else:
        query_ref_index_label_alignments = Index_label_alignments(blastn_alignments_lines)
        count = 1
        while count > 0:
            count = query_ref_index_label_alignments.pop_contained_alignments()
        one_alignments = query_ref_index_label_alignments.get_sorted_one_alignments()
        to_print_1 = ""
        to_print_2 = ""
        # 0	1	ERR9808518.1002029;mito_1;98.724;12071;37;111;3;12027;87801;99800;0.0;21329
        # 1	1	ERR9808518.1002029;mito_1;99.269;7384;10;41;5672;13041;339447;332094;0.0;13296
        q_id = one_alignments[0][2].split(";")[0]
        s_id = one_alignments[0][2].split(";")[1]
        num_align = len(one_alignments)
        strand_list = [ "+" for i in range(num_align) ]
        percent_ident_list = [ 100 for i in range(num_align) ]
        qs_list = [1] * num_align
        qe_list = [1] * num_align
        ss_list = [1] * num_align
        se_list = [1] * num_align
        part_len_list = [1] * num_align
        AO_len_list = [1] * num_align
        copy_info_list = [""] * num_align # check for multiple mapping, disabled
        for i in range(num_align):
            align_info_fields = one_alignments[i][2].split(";")
            percent_ident_list[i] = float(align_info_fields[2])
            qs_list[i] = int(align_info_fields[6])
            qe_list[i] = int(align_info_fields[7])
            ss_list[i] = int(align_info_fields[8])
            se_list[i] = int(align_info_fields[9])
            part_len_list[i] = qe_list[i] - qs_list[i] + 1
            if ss_list[i] > se_list[i]:
                strand_list[i] = "-"
        for i in range(num_align):
            if (i+1) != num_align:
                AO_len_list[i] = qe_list[i] - qs_list[i+1] + 1
            else:
                AO_len_list[i] = "NA"
        percent_total = float((sum(part_len_list) - sum(AO_len_list[:-1])) / q_len * 100)
        # gather results
        aln_type = "aln_type=" + str(num_align) + ";"
        for i in range(num_align):
            if AO_len_list[i] == "NA":
                olp_type = "NA"
            elif type(AO_len_list[i]) == int:
                if AO_len_list[i] == 0:
                    olp_type = "ref"
                elif AO_len_list[i] > 0:
                    olp_type = "rep"
                elif AO_len_list[i] < 0:
                    olp_type = "ins"
            if i != (num_align - 1):
                aln_type += olp_type + ","
            else:
                aln_type += olp_type
        to_print_1 = "\t".join([q_id, s_id, str(q_len), aln_type, str(percent_total)])
        # check if each part is dual mapping, disabled
        for i in range(num_align):
            copy_info_list[i] = "cn=1;c1=100.0,1,1"
            # gather results
            to_print_2 += "\taln=" + str(i+1) + ";len=" + str(part_len_list[i]) + ";olp=" + str(AO_len_list[i]) + ";idt=" + str(percent_ident_list[i]) + ";strand=" + strand_list[i] + ";qs=" + str(qs_list[i]) + ";qe=" + str(qe_list[i]) + ";ss=" + str(ss_list[i]) + ";se=" + str(se_list[i]) + ";" + copy_info_list[i]
        with open(tmp_blastn_dir + "/" + read_record.id + "_sorted_blastn_alignments.txt", "wt") as fout:
            print(to_print_1 + to_print_2, file=fout)
    command_1 = "rm "+ tmp_blastn_dir + "/" + read_record.id + ".fasta"
    command_2 = "cat "+ tmp_blastn_dir + "/" + read_record.id + "_sorted_blastn_alignments.txt >> " + tmp_blastn_dir + "/" + results_filename
    command_3 = "rm "+ tmp_blastn_dir + "/" + read_record.id + "_sorted_blastn_alignments.txt"
    commands = command_1 + ";" + command_2 + ";" + command_3
    hfbase.get_cli_output_lines(commands, side_effect = True)
    return


def run_multi_threads_blastn(sample_index, genome, run_info, reads_filename, genome_absolute_path, results_filename, soft_paths_dict, threads):
    if "/" in results_filename:
        print("Error: results_filename contains '/'", file=sys.stderr)
        sys.exit(1)
    tmp_root = get_tmp_root()
    tmp_blastn_dir = os.path.join(tmp_root, sample_index, genome, run_info, "tmp_blastn_results")
    if not os.path.exists(tmp_blastn_dir):
        os.makedirs(tmp_blastn_dir)
    if os.path.exists(tmp_blastn_dir + "/" + results_filename):
        os.remove(tmp_blastn_dir + "/" + results_filename)
    command_1 = soft_paths_dict.get("seqkit") + " fq2fa " + reads_filename + " > reads.fasta"
    hfbase.get_cli_output_lines(command_1, side_effect = True)
    hfreads.replace_reads_id("reads.fasta", "new_reads.fasta")
    read_records = SeqIO.parse("new_reads.fasta", "fasta")
    if os.path.exists(results_filename):
        os.remove(results_filename)
    with cf.ThreadPoolExecutor(int(threads)) as tex:
        futures = [tex.submit(run_blastn_sorter_single, read_record, genome_absolute_path, results_filename, soft_paths_dict, tmp_blastn_dir) for read_record in read_records]
        results = [future.result() for future in cf.as_completed(futures)]
    command_1 = "mv " + tmp_blastn_dir + "/" + results_filename + " ."
    command_2 = "rm -rf " + os.path.join(tmp_root, sample_index)
    commands = command_1 + ";" + command_2
    hfbase.get_cli_output_lines(commands, side_effect = True)


def get_type_and_subtype(blastn_info_file, default_num_types, out_dir="read_group_files"):
    blastn_results = hfbase.get_file_lines(blastn_info_file)
    default_num_types = 5 # type_1,2,3,4,5 and other
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    for line in blastn_results:
        fields = line.split("\t")
        type_info = fields[3]
        aln_type, aln_subtype = type_info.split(";") # aln_type=4;rep,rep,rep,NA
        num_align = int(aln_type.split("=")[1])
        if num_align > default_num_types:
            with open(out_dir + "/other_blastn_results.txt", "at") as fout:
                print(line, file=fout)
        else:
            subtype = "_".join(aln_subtype.split(","))
            out_file = out_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_blastn_results.txt"
            with open(out_file, "at") as fout:
                print(line, file=fout)


def check_FL_and_multi(blastn_info_file, default_num_types, out_dir="FL_read_group_files", id_dir="IDs", report_dir="Reports"):
    blastn_df = pd.read_table(blastn_info_file, header=None)
    # get the type and subtype
    type_info_list = list(set(blastn_df.iloc[:,3]))
    if len(type_info_list) == 1:
        type_info = type_info_list[0] # aln_type=1;NA
        aln_type, aln_subtype = type_info.split(";") # aln_type=4;rep,rep,rep,NA
        num_align = int(aln_type.split("=")[1])
        subtype = "_".join(aln_subtype.split(",")) # rep_rep_rep_NA
    else:
        print("Error: more than one type info in blastn file")
        for type_info in type_info_list:
            print(type_info)
        sys.exit(1)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if not os.path.exists(id_dir):
        os.mkdir(id_dir)
    if not os.path.exists(report_dir):
        os.mkdir(report_dir)
    if num_align > default_num_types: # only check for type_1,2,3,4,5
        return
    count_FL = 0
    count_FL_multi = 0
    count_partial = 0
    for i in range(len(blastn_df)):
        read_id = blastn_df.iloc[i,0]
        percent_total = blastn_df.iloc[i,4]
        align_info_list_of_dict = list()
        for j in range(0,num_align):
            tmp_aln_info_list = blastn_df.iloc[i,5+j].split(";") # ['aln=1', 'len=16011', 'olp=NA', 'idt=99.975', 'strand=+', 'qs=1', 'qe=16011', 'ss=73154', 'se=89162', 'cn=1', 'c1=100.0,73154,89162']
            tmp_od = OrderedDict()
            for tmp_aln_info in tmp_aln_info_list:
                tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
                tmp_od[tmp_aln_info_key] = tmp_aln_info_value
            align_info_list_of_dict.append(tmp_od)
        # check if the read is full length
        if percent_total >= 98: # get IDs of full length reads (98%) and partial reads
            with open(id_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_FL_ids.txt", "at") as fout:
                print(read_id, file=fout)
                count_FL += 1
            multi_flag = False
            for j in range(0,num_align):
                if int(align_info_list_of_dict[j]["cn"]) > 1:
                    multi_flag = True
                    break
            if multi_flag:
                with open(id_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_FL_multi_ids.txt", "at") as fout:
                    print(read_id, file=fout)
                    count_FL_multi += 1
        else:
            with open(id_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_partial_ids.txt", "at") as fout:
                print(read_id, file=fout)
                count_partial += 1
    # report the number of reads in each group
    with open(report_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_report.txt", "wt") as fout:
        print("FL", count_FL, sep="\t", file=fout)
        print("multiple mapping", count_FL_multi, sep="\t", file=fout)
        print("partial", count_partial, sep="\t", file=fout)
    # get the blastn info of full length reads from the original blastn file
    blastn_results = hfbase.get_file_lines(blastn_info_file)
    with open(out_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_FL_blastn_results.txt", "wt") as fout:
        for i in range(len(blastn_results)):
            if float(blastn_results[i].split("\t")[4]) > 98:
                print(blastn_results[i], file=fout)


def get_next_groups(groups):
    next_groups = list()
    tmp_groups_ins = [["ins"] + group for group in groups]
    tmp_groups_ref = [["ref"] + group for group in groups]
    tmp_groups_rep = [["rep"] + group for group in groups]
    next_groups = tmp_groups_ins + tmp_groups_ref + tmp_groups_rep
    return next_groups


def match_se1_ss2(old_index, num_align, subtype, SE1, SS2, blastn_df, ID_subgroup_dir="ID_subgroup"):
    subgroup_df = pd.DataFrame(columns=["read_id", "olp_1", "strand_1", "ss_1", "se_1", "cn_1", "olp_2", "strand_2", "ss_2", "se_2", "cn_2"])
    subgroup_count = 0
    for i in range(len(blastn_df)):
        # parse the info
        read_id = blastn_df.iloc[i,0]
        align_info_list_of_dict = list()
        for j in range(0, num_align):
            tmp_aln_info_list = blastn_df.iloc[i,5+j].split(";") # ['aln=1', 'len=16011', 'olp=NA', 'idt=99.975', 'strand=+', 'qs=1', 'qe=16011', 'ss=73154', 'se=89162', 'cn=1', 'c1=100.0,73154,89162']
            tmp_od = OrderedDict()
            for tmp_aln_info in tmp_aln_info_list:
                tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
                tmp_od[tmp_aln_info_key] = tmp_aln_info_value
            align_info_list_of_dict.append(tmp_od)
        # extract the info
        aln_len_1 = int(align_info_list_of_dict[0]["len"])
        olp_1 = align_info_list_of_dict[0]["olp"]
        idt_1 = float(align_info_list_of_dict[0]["idt"])
        strand_1 = align_info_list_of_dict[0]["strand"]
        ss_1 = int(align_info_list_of_dict[0]["ss"])
        se_1 = int(align_info_list_of_dict[0]["se"])
        cn_1 = int(align_info_list_of_dict[0]["cn"])
        aln_len_2 = int(align_info_list_of_dict[1]["len"])
        olp_2 = align_info_list_of_dict[1]["olp"] # NA
        idt_2 = float(align_info_list_of_dict[1]["idt"])
        strand_2 = align_info_list_of_dict[1]["strand"]
        ss_2 = int(align_info_list_of_dict[1]["ss"])
        se_2 = int(align_info_list_of_dict[1]["se"])
        cn_2 = int(align_info_list_of_dict[1]["cn"])
        # generate a new dataframe with columns: read_id, olp_1, ss_1, se_1, cn_1, olp_2, ss_2, se_2, cn_2 
        if se_1 == SE1 and ss_2 == SS2: # update the dataframe when se_1 == SE1 and ss_2 == SS2
            subgroup_df.loc[subgroup_count] = [read_id, olp_1, strand_1, ss_1, se_1, cn_1, olp_2, strand_2, ss_2, se_2, cn_2]
            subgroup_count += 1
    # analyze the dataframe, and write the info to a line in a combined summary dataframe
    fout_1 = open(ID_subgroup_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_subgroup_" + str(old_index) + "_FL_ids.txt", "wt")
    fout_2 = open(ID_subgroup_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_subgroup_" + str(old_index) + "_FL_multi_ids.txt", "wt")
    min_olp_1 = int(subgroup_df.iloc[0,1])
    mid_olp_1 = float(subgroup_df.iloc[0,1])
    max_olp_1 = int(subgroup_df.iloc[0,1])
    subgroup_count = len(subgroup_df)
    subgroup_multi_count = 0
    for i in range(len(subgroup_df)):
        # [read_id, olp_1, strand_1, ss_1, se_1, cn_1, olp_2, strand_2, ss_2, se_2, cn_2]
        read_id = subgroup_df.iloc[i,0]
        olp_1 = int(subgroup_df.iloc[i,1])
        strand_1 = subgroup_df.iloc[i,2]
        ss_1 = subgroup_df.iloc[i,3]
        se_1 = subgroup_df.iloc[i,4]
        cn_1 = subgroup_df.iloc[i,5]
        olp_2 = subgroup_df.iloc[i,6]
        strand_2 = subgroup_df.iloc[i,7]
        ss_2 = subgroup_df.iloc[i,8]
        se_2 = subgroup_df.iloc[i,9]
        cn_2 = subgroup_df.iloc[i,10]
        if olp_1 < min_olp_1:
            min_olp_1 = olp_1
        if olp_1 > max_olp_1:
            max_olp_1 = olp_1
        print(read_id, file=fout_1)
        if cn_1 > 1 or cn_2 > 1:
            subgroup_multi_count += 1
            print(read_id, file=fout_2)
    fout_1.close()
    fout_2.close()
    mid_olp_1 = (min_olp_1 + max_olp_1) / 2
    strand_1_most = subgroup_df["strand_1"].mode()[0]
    strand_2_most = subgroup_df["strand_2"].mode()[0]
    return min_olp_1, mid_olp_1, max_olp_1, strand_1_most, strand_2_most, subgroup_count, subgroup_multi_count


def match_se1_ss2_se2_ss3(old_index, num_align, subtype, SE1, SS2, SE2, SS3, blastn_df, ID_subgroup_dir="ID_subgroup"):
    subgroup_df = pd.DataFrame(columns=["read_id", "olp_1", "strand_1", "ss_1", "se_1", "cn_1", "olp_2", "strand_2", "ss_2", "se_2", "cn_2", "olp_3", "strand_3", "ss_3", "se_3", "cn_3"])
    subgroup_count = 0
    for i in range(len(blastn_df)):
        # parse the info
        read_id = blastn_df.iloc[i,0]
        align_info_list_of_dict = list()
        for j in range(0, num_align):
            tmp_aln_info_list = blastn_df.iloc[i,5+j].split(";") # ['aln=1', 'len=16011', 'olp=NA', 'idt=99.975', 'strand=+', 'qs=1', 'qe=16011', 'ss=73154', 'se=89162', 'cn=1', 'c1=100.0,73154,89162']
            tmp_od = OrderedDict()
            for tmp_aln_info in tmp_aln_info_list:
                tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
                tmp_od[tmp_aln_info_key] = tmp_aln_info_value
            align_info_list_of_dict.append(tmp_od)
        # extract the info
        aln_len_1 = int(align_info_list_of_dict[0]["len"])
        olp_1 = align_info_list_of_dict[0]["olp"]
        idt_1 = float(align_info_list_of_dict[0]["idt"])
        strand_1 = align_info_list_of_dict[0]["strand"]
        ss_1 = int(align_info_list_of_dict[0]["ss"])
        se_1 = int(align_info_list_of_dict[0]["se"])
        cn_1 = int(align_info_list_of_dict[0]["cn"])
        aln_len_2 = int(align_info_list_of_dict[1]["len"])
        olp_2 = align_info_list_of_dict[1]["olp"] # NA
        idt_2 = float(align_info_list_of_dict[1]["idt"])
        strand_2 = align_info_list_of_dict[1]["strand"]
        ss_2 = int(align_info_list_of_dict[1]["ss"])
        se_2 = int(align_info_list_of_dict[1]["se"])
        cn_2 = int(align_info_list_of_dict[1]["cn"])
        aln_len_3 = int(align_info_list_of_dict[2]["len"])
        olp_3 = align_info_list_of_dict[2]["olp"]
        idt_3 = float(align_info_list_of_dict[2]["idt"])
        strand_3 = align_info_list_of_dict[2]["strand"]
        ss_3 = int(align_info_list_of_dict[2]["ss"])
        se_3 = int(align_info_list_of_dict[2]["se"])
        cn_3 = int(align_info_list_of_dict[2]["cn"])
        if se_1 == SE1 and ss_2 == SS2 and se_2 == SE2 and ss_3 == SS3: # update the dataframe when se_1 == SE1 and ss_2 == SS2 and se_2 == SE2 and ss_3 == SS3
            subgroup_df.loc[subgroup_count] = [read_id, olp_1, strand_1, ss_1, se_1, cn_1, olp_2, strand_2, ss_2, se_2, cn_2, olp_3, strand_3, ss_3, se_3, cn_3]
            subgroup_count += 1
    # analyze the dataframe, and write the info to a line in a combined summary dataframe
    fout_1 = open(ID_subgroup_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_subgroup_" + str(old_index) + "_FL_ids.txt", "wt")
    fout_2 = open(ID_subgroup_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_subgroup_" + str(old_index) + "_FL_multi_ids.txt", "wt")
    min_olp_1 = int(subgroup_df.iloc[0,1])
    mid_olp_1 = float(subgroup_df.iloc[0,1])
    max_olp_1 = int(subgroup_df.iloc[0,1])
    min_olp_2 = int(subgroup_df.iloc[0,5])
    mid_olp_2 = float(subgroup_df.iloc[0,5])
    max_olp_2 = int(subgroup_df.iloc[0,5])
    # ["read_id", "olp_1", "ss_1", "se_1", "cn_1", "olp_2", "ss_2", "se_2", "cn_2", "olp_3", "ss_3", "se_3", "cn_3"]
    subgroup_count = len(subgroup_df)
    subgroup_multi_count = 0
    for i in range(len(subgroup_df)):
        read_id = subgroup_df.iloc[i,0]
        olp_1 = int(subgroup_df.iloc[i,1])
        strand_1 = subgroup_df.iloc[i,2]
        ss_1 = subgroup_df.iloc[i,3]
        se_1 = subgroup_df.iloc[i,4]
        cn_1 = subgroup_df.iloc[i,5]
        olp_2 = int(subgroup_df.iloc[i,6])
        strand_2 = subgroup_df.iloc[i,7]
        ss_2 = subgroup_df.iloc[i,8]
        se_2 = subgroup_df.iloc[i,9]
        cn_2 = subgroup_df.iloc[i,10]
        olp_3 = subgroup_df.iloc[i,11]
        strand_3 = subgroup_df.iloc[i,12]
        ss_3 = subgroup_df.iloc[i,13]
        se_3 = subgroup_df.iloc[i,14]
        cn_3 = subgroup_df.iloc[i,15]
        if olp_1 < min_olp_1:
            min_olp_1 = olp_1
        if olp_1 > max_olp_1:
            max_olp_1 = olp_1
        if olp_2 < min_olp_2:
            min_olp_2 = olp_2
        if olp_2 > max_olp_2:
            max_olp_2 = olp_2
        print(read_id, file=fout_1)
        if cn_1 > 1 or cn_2 > 1 or cn_3 > 1:
            subgroup_multi_count += 1
            print(read_id, file=fout_2)
    fout_1.close()
    fout_2.close()
    mid_olp_1 = (min_olp_1 + max_olp_1) / 2
    mid_olp_2 = (min_olp_2 + max_olp_2) / 2
    strand_1_most = subgroup_df["strand_1"].mode()[0]
    strand_2_most = subgroup_df["strand_2"].mode()[0]
    strand_3_most = subgroup_df["strand_3"].mode()[0]
    return min_olp_1, mid_olp_1, max_olp_1, min_olp_2, mid_olp_2, max_olp_2, strand_1_most, strand_2_most, strand_3_most, subgroup_count, subgroup_multi_count


def match_se1_ss2_se2_ss3_se3_ss4(old_index, num_align, subtype, SE1, SS2, SE2, SS3, SE3, SS4, blastn_df, ID_subgroup_dir="ID_subgroup"):
    subgroup_df = pd.DataFrame(columns=["read_id", "olp_1", "strand_1", "ss_1", "se_1", "cn_1", "olp_2", "strand_2", "ss_2", "se_2", "cn_2", "olp_3", "strand_3", "ss_3", "se_3", "cn_3", "olp_4", "strand_4", "ss_4", "se_4", "cn_4"])
    subgroup_count = 0
    for i in range(len(blastn_df)):
        # parse the info
        read_id = blastn_df.iloc[i,0]
        align_info_list_of_dict = list()
        for j in range(0, num_align):
            tmp_aln_info_list = blastn_df.iloc[i,5+j].split(";") # ['aln=1', 'len=16011', 'olp=NA', 'idt=99.975', 'strand=+', 'qs=1', 'qe=16011', 'ss=73154', 'se=89162', 'cn=1', 'c1=100.0,73154,89162']
            tmp_od = OrderedDict()
            for tmp_aln_info in tmp_aln_info_list:
                tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
                tmp_od[tmp_aln_info_key] = tmp_aln_info_value
            align_info_list_of_dict.append(tmp_od)
        # extract the info
        aln_len_1 = int(align_info_list_of_dict[0]["len"])
        olp_1 = align_info_list_of_dict[0]["olp"]
        idt_1 = float(align_info_list_of_dict[0]["idt"])
        strand_1 = align_info_list_of_dict[0]["strand"]
        ss_1 = int(align_info_list_of_dict[0]["ss"])
        se_1 = int(align_info_list_of_dict[0]["se"])
        cn_1 = int(align_info_list_of_dict[0]["cn"])
        aln_len_2 = int(align_info_list_of_dict[1]["len"])
        olp_2 = align_info_list_of_dict[1]["olp"] # NA
        idt_2 = float(align_info_list_of_dict[1]["idt"])
        strand_2 = align_info_list_of_dict[1]["strand"]
        ss_2 = int(align_info_list_of_dict[1]["ss"])
        se_2 = int(align_info_list_of_dict[1]["se"])
        cn_2 = int(align_info_list_of_dict[1]["cn"])
        aln_len_3 = int(align_info_list_of_dict[2]["len"])
        olp_3 = align_info_list_of_dict[2]["olp"]
        idt_3 = float(align_info_list_of_dict[2]["idt"])
        strand_3 = align_info_list_of_dict[2]["strand"]
        ss_3 = int(align_info_list_of_dict[2]["ss"])
        se_3 = int(align_info_list_of_dict[2]["se"])
        cn_3 = int(align_info_list_of_dict[2]["cn"])
        aln_len_4 = int(align_info_list_of_dict[3]["len"])
        olp_4 = align_info_list_of_dict[3]["olp"]
        idt_4 = float(align_info_list_of_dict[3]["idt"])
        strand_4 = align_info_list_of_dict[3]["strand"]
        ss_4 = int(align_info_list_of_dict[3]["ss"])
        se_4 = int(align_info_list_of_dict[3]["se"])
        cn_4 = int(align_info_list_of_dict[3]["cn"])
        if se_1 == SE1 and ss_2 == SS2 and se_2 == SE2 and ss_3 == SS3 and se_3 == SE3 and ss_4 == SS4: # update the dataframe when se_1 == SE1 and ss_2 == SS2 and se_2 == SE2 and ss_3 == SS3 and se_3 == SE3 and ss_4 == SS4
            subgroup_df.loc[subgroup_count] = [read_id, olp_1, strand_1, ss_1, se_1, cn_1, olp_2, strand_2, ss_2, se_2, cn_2, olp_3, strand_3, ss_3, se_3, cn_3, olp_4, strand_4, ss_4, se_4, cn_4]
            subgroup_count += 1
    # analyze the dataframe, and write the info to a line in a combined summary dataframe
    fout_1 = open(ID_subgroup_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_subgroup_" + str(old_index) + "_FL_ids.txt", "wt")
    fout_2 = open(ID_subgroup_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_subgroup_" + str(old_index) + "_FL_multi_ids.txt", "wt")
    min_olp_1 = int(subgroup_df.iloc[0,1])
    mid_olp_1 = float(subgroup_df.iloc[0,1])
    max_olp_1 = int(subgroup_df.iloc[0,1])
    min_olp_2 = int(subgroup_df.iloc[0,5])
    mid_olp_2 = float(subgroup_df.iloc[0,5])
    max_olp_2 = int(subgroup_df.iloc[0,5])
    min_olp_3 = int(subgroup_df.iloc[0,9])
    mid_olp_3 = float(subgroup_df.iloc[0,9])
    max_olp_3 = int(subgroup_df.iloc[0,9])
    # ["read_id", "olp_1", "ss_1", "se_1", "cn_1", "olp_2", "ss_2", "se_2", "cn_2", "olp_3", "ss_3", "se_3", "cn_3", "olp_4", "ss_4", "se_4", "cn_4"]
    subgroup_count = len(subgroup_df)
    subgroup_multi_count = 0
    for i in range(len(subgroup_df)):
        read_id = subgroup_df.iloc[i,0]
        olp_1 = int(subgroup_df.iloc[i,1])
        strand_1 = subgroup_df.iloc[i,2]
        ss_1 = subgroup_df.iloc[i,3]
        se_1 = subgroup_df.iloc[i,4]
        cn_1 = subgroup_df.iloc[i,5]
        olp_2 = int(subgroup_df.iloc[i,6])
        strand_2 = subgroup_df.iloc[i,7]
        ss_2 = subgroup_df.iloc[i,8]
        se_2 = subgroup_df.iloc[i,9]
        cn_2 = subgroup_df.iloc[i,10]
        olp_3 = int(subgroup_df.iloc[i,11])
        strand_3 = subgroup_df.iloc[i,12]
        ss_3 = subgroup_df.iloc[i,13]
        se_3 = subgroup_df.iloc[i,14]
        cn_3 = subgroup_df.iloc[i,15]
        olp_4 = subgroup_df.iloc[i,16]
        strand_4 = subgroup_df.iloc[i,17]
        ss_4 = subgroup_df.iloc[i,18]
        se_4 = subgroup_df.iloc[i,19]
        cn_4 = subgroup_df.iloc[i,20]
        if olp_1 < min_olp_1:
            min_olp_1 = olp_1
        if olp_1 > max_olp_1:
            max_olp_1 = olp_1
        if olp_2 < min_olp_2:
            min_olp_2 = olp_2
        if olp_2 > max_olp_2:
            max_olp_2 = olp_2
        if olp_3 < min_olp_3:
            min_olp_3 = olp_3
        if olp_3 > max_olp_3:
            max_olp_3 = olp_3
        print(read_id, file=fout_1)
        if cn_1 > 1 or cn_2 > 1 or cn_3 > 1 or cn_4 > 1:
            subgroup_multi_count += 1
            print(read_id, file=fout_2)
    fout_1.close()
    fout_2.close()
    mid_olp_1 = (min_olp_1 + max_olp_1) / 2
    mid_olp_2 = (min_olp_2 + max_olp_2) / 2
    mid_olp_3 = (min_olp_3 + max_olp_3) / 2
    strand_1_most = subgroup_df["strand_1"].mode()[0]
    strand_2_most = subgroup_df["strand_2"].mode()[0]
    strand_3_most = subgroup_df["strand_3"].mode()[0]
    strand_4_most = subgroup_df["strand_4"].mode()[0]
    return min_olp_1, mid_olp_1, max_olp_1, min_olp_2, mid_olp_2, max_olp_2, min_olp_3, mid_olp_3, max_olp_3, strand_1_most, strand_2_most, strand_3_most, strand_4_most, subgroup_count, subgroup_multi_count


def match_se1_ss2_se2_ss3_se3_ss4_se4_ss5(old_index, num_align, subtype, SE1, SS2, SE2, SS3, SE3, SS4, SE4, SS5, blastn_df, ID_subgroup_dir="ID_subgroup"):
    subgroup_df = pd.DataFrame(columns=["read_id", "olp_1", "strand_1", "ss_1", "se_1", "cn_1", "olp_2", "strand_2", "ss_2", "se_2", "cn_2", "olp_3", "strand_3", "ss_3", "se_3", "cn_3", "olp_4", "strand_4", "ss_4", "se_4", "cn_4", "olp_5", "strand_5", "ss_5", "se_5", "cn_5"])
    subgroup_count = 0
    for i in range(len(blastn_df)):
        # parse the info
        read_id = blastn_df.iloc[i,0]
        align_info_list_of_dict = list()
        for j in range(0, num_align):
            tmp_aln_info_list = blastn_df.iloc[i,5+j].split(";") # ['aln=1', 'len=16011', 'olp=NA', 'idt=99.975', 'strand=+', 'qs=1', 'qe=16011', 'ss=73154', 'se=89162', 'cn=1', 'c1=100.0,73154,89162']
            tmp_od = OrderedDict()
            for tmp_aln_info in tmp_aln_info_list:
                tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
                tmp_od[tmp_aln_info_key] = tmp_aln_info_value
            align_info_list_of_dict.append(tmp_od)
        # extract the info
        aln_len_1 = int(align_info_list_of_dict[0]["len"])
        olp_1 = align_info_list_of_dict[0]["olp"]
        idt_1 = float(align_info_list_of_dict[0]["idt"])
        strand_1 = align_info_list_of_dict[0]["strand"]
        ss_1 = int(align_info_list_of_dict[0]["ss"])
        se_1 = int(align_info_list_of_dict[0]["se"])
        cn_1 = int(align_info_list_of_dict[0]["cn"])
        aln_len_2 = int(align_info_list_of_dict[1]["len"])
        olp_2 = align_info_list_of_dict[1]["olp"] # NA
        idt_2 = float(align_info_list_of_dict[1]["idt"])
        strand_2 = align_info_list_of_dict[1]["strand"]
        ss_2 = int(align_info_list_of_dict[1]["ss"])
        se_2 = int(align_info_list_of_dict[1]["se"])
        cn_2 = int(align_info_list_of_dict[1]["cn"])
        aln_len_3 = int(align_info_list_of_dict[2]["len"])
        olp_3 = align_info_list_of_dict[2]["olp"]
        idt_3 = float(align_info_list_of_dict[2]["idt"])
        strand_3 = align_info_list_of_dict[2]["strand"]
        ss_3 = int(align_info_list_of_dict[2]["ss"])
        se_3 = int(align_info_list_of_dict[2]["se"])
        cn_3 = int(align_info_list_of_dict[2]["cn"])
        aln_len_4 = int(align_info_list_of_dict[3]["len"])
        olp_4 = align_info_list_of_dict[3]["olp"]
        idt_4 = float(align_info_list_of_dict[3]["idt"])
        strand_4 = align_info_list_of_dict[3]["strand"]
        ss_4 = int(align_info_list_of_dict[3]["ss"])
        se_4 = int(align_info_list_of_dict[3]["se"])
        cn_4 = int(align_info_list_of_dict[3]["cn"])
        aln_len_5 = int(align_info_list_of_dict[4]["len"])
        olp_5 = align_info_list_of_dict[4]["olp"]
        idt_5 = float(align_info_list_of_dict[4]["idt"])
        strand_5 = align_info_list_of_dict[4]["strand"]
        ss_5 = int(align_info_list_of_dict[4]["ss"])
        se_5 = int(align_info_list_of_dict[4]["se"])
        cn_5 = int(align_info_list_of_dict[4]["cn"])
        if se_1 == SE1 and ss_2 == SS2 and se_2 == SE2 and ss_3 == SS3 and se_3 == SE3 and ss_4 == SS4 and se_4 == SE4 and ss_5 == SS5: # update the dataframe when se_1 == SE1 and ss_2 == SS2 and se_2 == SE2 and ss_3 == SS3 and se_3 == SE3 and ss_4 == SS4 and se_4 == SE4 and ss_5 == SS5
            subgroup_df.loc[subgroup_count] = [read_id, olp_1, strand_1, ss_1, se_1, cn_1, olp_2, strand_2, ss_2, se_2, cn_2, olp_3, strand_3, ss_3, se_3, cn_3, olp_4, strand_4, ss_4, se_4, cn_4, olp_5, strand_5, ss_5, se_5, cn_5]
            subgroup_count += 1
    # analyze the dataframe, and write the info to a line in a combined summary dataframe
    fout_1 = open(ID_subgroup_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_subgroup_" + str(old_index) + "_FL_ids.txt", "wt")
    fout_2 = open(ID_subgroup_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_subgroup_" + str(old_index) + "_FL_multi_ids.txt", "wt")
    min_olp_1 = int(subgroup_df.iloc[0,1])
    mid_olp_1 = float(subgroup_df.iloc[0,1])
    max_olp_1 = int(subgroup_df.iloc[0,1])
    min_olp_2 = int(subgroup_df.iloc[0,5])
    mid_olp_2 = float(subgroup_df.iloc[0,5])
    max_olp_2 = int(subgroup_df.iloc[0,5])
    min_olp_3 = int(subgroup_df.iloc[0,9])
    mid_olp_3 = float(subgroup_df.iloc[0,9])
    max_olp_3 = int(subgroup_df.iloc[0,9])
    min_olp_4 = int(subgroup_df.iloc[0,13])
    mid_olp_4 = float(subgroup_df.iloc[0,13])
    max_olp_4 = int(subgroup_df.iloc[0,13])
    # ["read_id", "olp_1", "ss_1", "se_1", "cn_1", "olp_2", "ss_2", "se_2", "cn_2", "olp_3", "ss_3", "se_3", "cn_3", "olp_4", "ss_4", "se_4", "cn_4", "olp_5", "ss_5", "se_5", "cn_5"]
    subgroup_count = len(subgroup_df)
    subgroup_multi_count = 0
    for i in range(len(subgroup_df)):
        read_id = subgroup_df.iloc[i,0]
        olp_1 = int(subgroup_df.iloc[i,1])
        strand_1 = subgroup_df.iloc[i,2]
        ss_1 = subgroup_df.iloc[i,3]
        se_1 = subgroup_df.iloc[i,4]
        cn_1 = subgroup_df.iloc[i,5]
        olp_2 = int(subgroup_df.iloc[i,6])
        strand_2 = subgroup_df.iloc[i,7]
        ss_2 = subgroup_df.iloc[i,8]
        se_2 = subgroup_df.iloc[i,9]
        cn_2 = subgroup_df.iloc[i,10]
        olp_3 = int(subgroup_df.iloc[i,11])
        strand_3 = subgroup_df.iloc[i,12]
        ss_3 = subgroup_df.iloc[i,13]
        se_3 = subgroup_df.iloc[i,14]
        cn_3 = subgroup_df.iloc[i,15]
        olp_4 = int(subgroup_df.iloc[i,16])
        strand_4 = subgroup_df.iloc[i,17]
        ss_4 = subgroup_df.iloc[i,18]
        se_4 = subgroup_df.iloc[i,19]
        cn_4 = subgroup_df.iloc[i,20]
        olp_5 = subgroup_df.iloc[i,21]
        strand_5 = subgroup_df.iloc[i,22]
        ss_5 = subgroup_df.iloc[i,23]
        se_5 = subgroup_df.iloc[i,24]
        cn_5 = subgroup_df.iloc[i,25]
        if olp_1 < min_olp_1:
            min_olp_1 = olp_1
        if olp_1 > max_olp_1:
            max_olp_1 = olp_1
        if olp_2 < min_olp_2:
            min_olp_2 = olp_2
        if olp_2 > max_olp_2:
            max_olp_2 = olp_2
        if olp_3 < min_olp_3:
            min_olp_3 = olp_3
        if olp_3 > max_olp_3:
            max_olp_3 = olp_3
        if olp_4 < min_olp_4:
            min_olp_4 = olp_4
        if olp_4 > max_olp_4:
            max_olp_4 = olp_4
        print(read_id, file=fout_1)
        if cn_1 > 1 or cn_2 > 1 or cn_3 > 1 or cn_4 > 1 or cn_5 > 1:
            subgroup_multi_count += 1
            print(read_id, file=fout_2)
    fout_1.close()
    fout_2.close()
    mid_olp_1 = (min_olp_1 + max_olp_1) / 2
    mid_olp_2 = (min_olp_2 + max_olp_2) / 2
    mid_olp_3 = (min_olp_3 + max_olp_3) / 2
    mid_olp_4 = (min_olp_4 + max_olp_4) / 2
    strand_1_most = subgroup_df["strand_1"].mode()[0]
    strand_2_most = subgroup_df["strand_2"].mode()[0]
    strand_3_most = subgroup_df["strand_3"].mode()[0]
    strand_4_most = subgroup_df["strand_4"].mode()[0]
    strand_5_most = subgroup_df["strand_5"].mode()[0]
    return min_olp_1, mid_olp_1, max_olp_1, min_olp_2, mid_olp_2, max_olp_2, min_olp_3, mid_olp_3, max_olp_3, min_olp_4, mid_olp_4, max_olp_4, strand_1_most, strand_2_most, strand_3_most, strand_4_most, strand_5_most, subgroup_count, subgroup_multi_count


def get_subgroups(groups, num_align, input_dir="FL_read_group_files", out_dir="combined_excel", ID_subgroup_dir="ID_subgroup"):
    if not os.path.exists(ID_subgroup_dir):
        os.mkdir(ID_subgroup_dir)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    for group in groups:
        subtype = "_".join(group)
        blastn_info_file = input_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_FL_blastn_results.txt"
        if os.path.exists(blastn_info_file) and os.path.getsize(blastn_info_file) > 0:
            blastn_df = pd.read_table(blastn_info_file, header=None)
            if num_align == 2:
                combine_df = pd.DataFrame(columns=["old_index", "(se1, ss2)", "strand_str", "min_olp_1", "mid_olp_1", "max_olp_1", "subgroup_count", "subgroup_multi_count"])
                all_se1_ss2 = set()
                for i in range(len(blastn_df)):
                    align_info_list_of_dict = list()
                    for j in range(0,num_align):
                        tmp_aln_info_list = blastn_df.iloc[i,5+j].split(";") # ['aln=1', 'len=16011', 'olp=NA', 'idt=99.975', 'strand=+', 'qs=1', 'qe=16011', 'ss=73154', 'se=89162', 'cn=1', 'c1=100.0,73154,89162']
                        tmp_od = OrderedDict()
                        for tmp_aln_info in tmp_aln_info_list:
                            tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
                            tmp_od[tmp_aln_info_key] = tmp_aln_info_value
                        align_info_list_of_dict.append(tmp_od)
                    all_se1_ss2.add((int(align_info_list_of_dict[0]["se"]), int(align_info_list_of_dict[1]["ss"])))
                all_se1_ss2_list = list(all_se1_ss2)
                all_se1_ss2_list.sort(key=lambda x:x[0]) # sort by se1
                indexed_all_se1_ss2 = [ (i+1, se1_ss2) for i, se1_ss2 in enumerate(all_se1_ss2_list)]
                for item in indexed_all_se1_ss2: # (1, (89162, 73154))
                    old_index = item[0]
                    SE1 = item[1][0] # a int
                    SS2 = item[1][1]
                    min_olp_1, mid_olp_1, max_olp_1, strand_1_most, strand_2_most, subgroup_count, subgroup_multi_count = match_se1_ss2(old_index, num_align, subtype, SE1, SS2, blastn_df, ID_subgroup_dir="ID_subgroup")
                    # ["new_index", (se1, ss2), min_olp_1, mid_olp_1, max_olp_1, "subgroup_count", "subgroup_multi_count"]
                    # update the combined summary dataframe, pay attention to datatypes
                    strand_str = strand_1_most + "," + strand_2_most
                    combine_df.loc[old_index-1] = [old_index, item[1], strand_str, min_olp_1, mid_olp_1, max_olp_1, subgroup_count, subgroup_multi_count]
                combine_df.sort_values(by="mid_olp_1", ascending=False, inplace=True) # sort combine_df by mid_olp_1
                # write the combined summary dataframe to a excel file
                combine_df.to_excel(out_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_summary_anno.xlsx", index=False)
                combine_df.to_excel(out_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_summary_backup.xlsx", index=False)
            elif num_align == 3:
                # get all (se1, ss2, se2, ss3)
                combine_df = pd.DataFrame(columns=["old_index", "(se1, ss2, se2, ss3)", "strand_str", "min_olp_1", "mid_olp_1", "max_olp_1", "min_olp_2", "mid_olp_2", "max_olp_2", "subgroup_count", "subgroup_multi_count"])
                all_se1_ss2_se2_ss3 = set()
                for i in range(len(blastn_df)):
                    align_info_list_of_dict = list()
                    for j in range(0,num_align):
                        tmp_aln_info_list = blastn_df.iloc[i,5+j].split(";") # ['aln=1', 'len=16011', 'olp=NA', 'idt=99.975', 'strand=+', 'qs=1', 'qe=16011', 'ss=73154', 'se=89162', 'cn=1', 'c1=100.0,73154,89162']
                        tmp_od = OrderedDict()
                        for tmp_aln_info in tmp_aln_info_list:
                            tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
                            tmp_od[tmp_aln_info_key] = tmp_aln_info_value
                        align_info_list_of_dict.append(tmp_od)
                    all_se1_ss2_se2_ss3.add((int(align_info_list_of_dict[0]["se"]), int(align_info_list_of_dict[1]["ss"]), int(align_info_list_of_dict[1]["se"]), int(align_info_list_of_dict[2]["ss"])))
                all_se1_ss2_se2_ss3_list = list(all_se1_ss2_se2_ss3)
                all_se1_ss2_se2_ss3_list.sort(key=lambda x:x[0]) # sort by se1
                indexed_all_se1_ss2_se2_ss3 = [ (i+1, se1_ss2_se2_ss3) for i, se1_ss2_se2_ss3 in enumerate(all_se1_ss2_se2_ss3_list)]
                for item in indexed_all_se1_ss2_se2_ss3: # (1, (89162, 73154, 89162, 73154))
                    old_index = item[0]
                    SE1 = item[1][0] # a int
                    SS2 = item[1][1]
                    SE2 = item[1][2]
                    SS3 = item[1][3]
                    min_olp_1, mid_olp_1, max_olp_1, min_olp_2, mid_olp_2, max_olp_2, strand_1_most, strand_2_most, strand_3_most, subgroup_count, subgroup_multi_count = match_se1_ss2_se2_ss3(old_index, num_align, subtype, SE1, SS2, SE2, SS3, blastn_df, ID_subgroup_dir="ID_subgroup")
                    # update the combined summary dataframe, pay attention to datatypes
                    strand_str = strand_1_most + "," + strand_2_most + "," + strand_3_most
                    combine_df.loc[old_index-1] = [old_index, item[1], strand_str, min_olp_1, mid_olp_1, max_olp_1, min_olp_2, mid_olp_2, max_olp_2, subgroup_count, subgroup_multi_count]
                combine_df.sort_values(by=["mid_olp_1", "mid_olp_2"], ascending=False, inplace=True) # sort combine_df by mid_olp_1 and mid_olp_2
                # write the combined summary dataframe to a excel file
                combine_df.to_excel(out_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_summary_anno.xlsx", index=False)
                combine_df.to_excel(out_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_summary_backup.xlsx", index=False)
            elif num_align == 4:
                # get all (se1, ss2, se2, ss3, se3, ss4)
                combine_df = pd.DataFrame(columns=["old_index", "(se1, ss2, se2, ss3)", "strand_str", "min_olp_1", "mid_olp_1", "max_olp_1", "min_olp_2", "mid_olp_2", "max_olp_2", "min_olp_3", "mid_olp_3", "max_olp_3", "subgroup_count", "subgroup_multi_count"])
                all_se1_ss2_se2_ss3_se3_ss4 = set()
                for i in range(len(blastn_df)):
                    align_info_list_of_dict = list()
                    for j in range(0,num_align):
                        tmp_aln_info_list = blastn_df.iloc[i,5+j].split(";") # ['aln=1', 'len=16011', 'olp=NA', 'idt=99.975', 'strand=+', 'qs=1', 'qe=16011', 'ss=73154', 'se=89162', 'cn=1', 'c1=100.0,73154,89162']
                        tmp_od = OrderedDict()
                        for tmp_aln_info in tmp_aln_info_list:
                            tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
                            tmp_od[tmp_aln_info_key] = tmp_aln_info_value
                        align_info_list_of_dict.append(tmp_od)
                    all_se1_ss2_se2_ss3_se3_ss4.add((int(align_info_list_of_dict[0]["se"]), int(align_info_list_of_dict[1]["ss"]), int(align_info_list_of_dict[1]["se"]), int(align_info_list_of_dict[2]["ss"]), int(align_info_list_of_dict[2]["se"]), int(align_info_list_of_dict[3]["ss"])))
                all_se1_ss2_se2_ss3_se3_ss4_list = list(all_se1_ss2_se2_ss3_se3_ss4)
                all_se1_ss2_se2_ss3_se3_ss4_list.sort(key=lambda x:x[0]) # sort by se1
                indexed_all_se1_ss2_se2_ss3_se3_ss4 = [ (i+1, se1_ss2_se2_ss3_se3_ss4) for i, se1_ss2_se2_ss3_se3_ss4 in enumerate(all_se1_ss2_se2_ss3_se3_ss4_list)]
                for item in indexed_all_se1_ss2_se2_ss3_se3_ss4: # (1, (89162, 73154, 89162, 73154, 89162, 73154))
                    old_index = item[0]
                    SE1 = item[1][0] # a int
                    SS2 = item[1][1]
                    SE2 = item[1][2]
                    SS3 = item[1][3]
                    SE3 = item[1][4]
                    SS4 = item[1][5]
                    min_olp_1, mid_olp_1, max_olp_1, min_olp_2, mid_olp_2, max_olp_2, min_olp_3, mid_olp_3, max_olp_3, strand_1_most, strand_2_most,strand_3_most, strand_4_most, subgroup_count, subgroup_multi_count = match_se1_ss2_se2_ss3_se3_ss4(old_index, num_align, subtype, SE1, SS2, SE2, SS3, SE3, SS4, blastn_df, ID_subgroup_dir="ID_subgroup")
                    # update the combined summary dataframe, pay attention to datatypes
                    strand_str = strand_1_most + "," + strand_2_most + "," + strand_3_most + "," + strand_4_most
                    combine_df.loc[old_index-1] = [old_index, item[1], strand_str, min_olp_1, mid_olp_1, max_olp_1, min_olp_2, mid_olp_2, max_olp_2, min_olp_3, mid_olp_3, max_olp_3, subgroup_count, subgroup_multi_count]
                combine_df.sort_values(by=["mid_olp_1", "mid_olp_2", "mid_olp_3"], ascending=False, inplace=True) # sort combine_df by mid_olp_1 and mid_olp_2 and mid_olp_3
                # write the combined summary dataframe to a excel file
                combine_df.to_excel(out_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_summary_anno.xlsx", index=False)
                combine_df.to_excel(out_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_summary_backup.xlsx", index=False)
            elif num_align == 5:
                # get all (se1, ss2, se2, ss3, se3, ss4, se4, ss5)
                combine_df = pd.DataFrame(columns=["old_index", "(se1, ss2, se2, ss3)", "strand_str", "min_olp_1", "mid_olp_1", "max_olp_1", "min_olp_2", "mid_olp_2", "max_olp_2", "min_olp_3", "mid_olp_3", "max_olp_3", "min_olp_4", "mid_olp_4", "max_olp_4", "subgroup_count", "subgroup_multi_count"])
                all_se1_ss2_se2_ss3_se3_ss4_se4_ss5 = set()
                for i in range(len(blastn_df)):
                    align_info_list_of_dict = list()
                    for j in range(0,num_align):
                        tmp_aln_info_list = blastn_df.iloc[i,5+j].split(";") # ['aln=1', 'len=16011', 'olp=NA', 'idt=99.975', 'strand=+', 'qs=1', 'qe=16011', 'ss=73154', 'se=89162', 'cn=1', 'c1=100.0,73154,89162']
                        tmp_od = OrderedDict()
                        for tmp_aln_info in tmp_aln_info_list:
                            tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
                            tmp_od[tmp_aln_info_key] = tmp_aln_info_value
                        align_info_list_of_dict.append(tmp_od)
                    all_se1_ss2_se2_ss3_se3_ss4_se4_ss5.add((int(align_info_list_of_dict[0]["se"]), int(align_info_list_of_dict[1]["ss"]), int(align_info_list_of_dict[1]["se"]), int(align_info_list_of_dict[2]["ss"]), int(align_info_list_of_dict[2]["se"]), int(align_info_list_of_dict[3]["ss"]), int(align_info_list_of_dict[3]["se"]), int(align_info_list_of_dict[4]["ss"])))
                all_se1_ss2_se2_ss3_se3_ss4_se4_ss5_list = list(all_se1_ss2_se2_ss3_se3_ss4_se4_ss5)
                all_se1_ss2_se2_ss3_se3_ss4_se4_ss5_list.sort(key=lambda x:x[0]) # sort by se1
                indexed_all_se1_ss2_se2_ss3_se3_ss4_se4_ss5 = [ (i+1, se1_ss2_se2_ss3_se3_ss4_se4_ss5) for i, se1_ss2_se2_ss3_se3_ss4_se4_ss5 in enumerate(all_se1_ss2_se2_ss3_se3_ss4_se4_ss5_list)]
                for item in indexed_all_se1_ss2_se2_ss3_se3_ss4_se4_ss5: # (1, (89162, 73154, 89162, 73154, 89162, 73154, 89162, 73154))
                    old_index = item[0]
                    SE1 = item[1][0] # a int
                    SS2 = item[1][1]
                    SE2 = item[1][2]
                    SS3 = item[1][3]
                    SE3 = item[1][4]
                    SS4 = item[1][5]
                    SE4 = item[1][6]
                    SS5 = item[1][7]
                    min_olp_1, mid_olp_1, max_olp_1, min_olp_2, mid_olp_2, max_olp_2, min_olp_3, mid_olp_3, max_olp_3, min_olp_4, mid_olp_4, max_olp_4, strand_1_most, strand_2_most,strand_3_most, strand_4_most, strand_5_most, subgroup_count, subgroup_multi_count = match_se1_ss2_se2_ss3_se3_ss4_se4_ss5(old_index, num_align, subtype, SE1, SS2, SE2, SS3, SE3, SS4, SE4, SS5, blastn_df, ID_subgroup_dir="ID_subgroup")
                    # update the combined summary dataframe, pay attention to datatypes
                    strand_str = strand_1_most + "," + strand_2_most + "," + strand_3_most + "," + strand_4_most + "," + strand_5_most
                    combine_df.loc[old_index-1] = [old_index, item[1], strand_str, min_olp_1, mid_olp_1, max_olp_1, min_olp_2, mid_olp_2, max_olp_2, min_olp_3, mid_olp_3, max_olp_3, min_olp_4, mid_olp_4, max_olp_4, subgroup_count, subgroup_multi_count]
                combine_df.sort_values(by=["mid_olp_1", "mid_olp_2", "mid_olp_3", "mid_olp_4"], ascending=False, inplace=True) # sort combine_df by mid_olp_1 and mid_olp_2 and mid_olp_3 and mid_olp_4
                # write the combined summary dataframe to a excel file
                combine_df.to_excel(out_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_summary_anno.xlsx", index=False)
                combine_df.to_excel(out_dir + "/type_" + str(num_align) + "_subtype_" + subtype + "_summary_backup.xlsx", index=False)
    return


def summarize_blastn_results(blastn_result_file="all_sorted_blastn_alignments.txt"):
    # get read groups
    get_type_and_subtype(blastn_result_file, 5, out_dir="read_group_files") # split files by type and subtype
    blastn_info_files = ["read_group_files/" + file for file in os.listdir("read_group_files") if file.startswith("type_")]
    for blastn_info_file in blastn_info_files:
        check_FL_and_multi(blastn_info_file, 5, out_dir="FL_read_group_files", id_dir="IDs", report_dir="Reports")
    # summarize the results
    FL_blastn_files = ["FL_read_group_files/" + file for file in os.listdir("FL_read_group_files") if file.endswith("_blastn_results.txt")]
    partial_blastn_files = ["IDs/" + file for file in os.listdir("IDs") if file.endswith("_partial_ids.txt")]
    other_blastn_file = "read_group_files/other_blastn_results.txt"
    if len(FL_blastn_files) == 0:
        command_1 = "touch all_FL_report.txt"
    else:
        command_1 = "wc -l " + " ".join(FL_blastn_files) + " > all_FL_report.txt"
    if len(partial_blastn_files) == 0:
        command_2 = "touch all_rm_report.txt"
    else:
        if not os.path.exists(other_blastn_file):
            command_2 = "wc -l " + " ".join(partial_blastn_files) + " > all_rm_report.txt"
        else:
            command_2 = "wc -l IDs/*_partial_ids.txt read_group_files/other_blastn_results.txt > all_rm_report.txt"
    command_3 = "wc -l " + blastn_result_file + " >> all_FL_report.txt"
    commands = command_1 + ";" + command_2 + ";" + command_3
    hfbase.get_cli_output_lines(commands, side_effect = True)
    # get read groups in excel format
    groups_1 = [["NA"]]
    groups_2 = get_next_groups(groups_1) # [["ins", "NA"], ["ref", "NA"], ["rep", "NA"]]
    groups_3 = get_next_groups(groups_2)
    groups_4 = get_next_groups(groups_3)
    groups_5 = get_next_groups(groups_4)
    groups_list = [groups_1, groups_2, groups_3, groups_4, groups_5] # for easy manipulation
    # analyze each subtype, and group the reads by (se1, ss2), (se1, ss2, se2, ss3), (se1, ss2, se2, ss3, se3, ss4), (se1, ss2, se2, ss3, se3, ss4, se4, ss5)
    for num_align in range(2, 6): # 2, 3, 4, 5
        get_subgroups(groups_list[num_align-1], num_align, input_dir="FL_read_group_files", out_dir="combined_excel", ID_subgroup_dir="ID_subgroup") # 1，2，3，4
    if  not os.path.exists("backup_info"):
        os.mkdir("backup_info")
    command_1 = "mv FL_read_group_files backup_info; mv IDs backup_info; mv Reports backup_info; mv read_group_files backup_info; mv " + blastn_result_file + " backup_info; mv all_rm_report.txt backup_info; mv combined_excel/*_backup.xlsx backup_info"
    command_2 = "find backup_info/FL_read_group_files -size 0 -delete"
    command_3 = "find ID_subgroup -size 0 -delete"
    commands = command_1 + ";" + command_2 + ";" + command_3
    hfbase.get_cli_output_lines(commands, side_effect = True)
    return


def get_cov_reads(current_dir, IDs_dir, soft_paths_dict, genome_absolute_path):
    reads_fasta = current_dir + "/new_reads.fasta" # clean_mito.fasta, not fastq
    ref_len = len(SeqIO.read(genome_absolute_path, "fasta").seq)
    # get type_2_ref_files, type_2_rep_large_files
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
    FL_ids_files = [IDs_dir + "/" + file for file in os.listdir(IDs_dir) if file.endswith("_FL_ids.txt")]
    partial_ids_files = [IDs_dir + "/" + file for file in os.listdir(IDs_dir) if file.endswith("_partial_ids.txt")]
    if len(FL_ids_files) == 0:
        command_1 = "touch FL_ids.txt"
    else:
        command_1 = "cat " + " ".join(FL_ids_files) + " > FL_ids.txt"
    if len(partial_ids_files) == 0:
        command_2 = "touch partial_ids.txt"
    else:
        command_2 = "cat " + " ".join(partial_ids_files) + " > partial_ids.txt"
    command_3 = soft_paths_dict.get("seqkit") + " grep -f FL_ids.txt " + reads_fasta + " > FL.fasta"
    command_4 = soft_paths_dict.get("seqkit") + " grep -f partial_ids.txt " + reads_fasta + " > partial.fasta"
    command_5 = "cut -f1 " + current_dir + "/backup_info/FL_read_group_files/type_1_subtype_NA_FL_blastn_results.txt " + " ".join(type_2_ref_files) + " " + " ".join(type_2_rep_large_files) + " > variant_cov_ids.txt"
    command_6 = soft_paths_dict.get("seqkit") + " grep -f variant_cov_ids.txt " + reads_fasta + " > variant_cov.fasta"
    commands = command_1 + " ; " + command_2 + " ; " + command_3 + " ; " + command_4 + " ; " + command_5 + " ; " + command_6
    hfbase.get_cli_output_lines(commands, side_effect = True)
    return ref_len


def run_bcftools(read_record, ref_fasta, sample_platform, soft_paths_dict, tmp_bcftools_dir):
    with open(tmp_bcftools_dir + "/" + read_record.id + ".fasta", "wt") as fout:
        print(read_record.format("fasta"), file=fout)
    platform_dict = {
            "HiFi": "map-hifi",
            "CLR": "map-pb",
            "ONT": "map-ont",
            "ultra-long": "map-ont",
        }
    command_1 = soft_paths_dict.get("minimap2") + " -ax " + platform_dict.get(sample_platform) + " " + ref_fasta + " " + tmp_bcftools_dir + "/" + read_record.id + ".fasta"
    command_2 = "samtools view -Sb -F 0x100 -@ 1 -"
    command_3 = "samtools sort -@ 1 -o -"
    command_4 = soft_paths_dict.get("bcftools") + " mpileup --indels-2.0 -m 1 -Ou -f " + ref_fasta + " -"
    command_5 = soft_paths_dict.get("bcftools") + " call -mv -P 0.99 -Ov | grep -v '^#' | cut -f2,4,5"
    command_6 = "rm " + tmp_bcftools_dir + "/" + read_record.id + ".fasta"
    commands = command_1 + " | " + command_2 + " | " + command_3 + " | " + command_4 + " | " + command_5 + " ; " + command_6
    results = hfbase.get_cli_output_lines(commands, side_effect = False)
    if len(results) > 0:
        with open(tmp_bcftools_dir + "/all_bcftools_calls.txt", "at") as fout:
            for result in results:
                fout.write(read_record.id + "\t" + result + "\n") # do not use print here, may lose "\n"
    return


def run_multi_threads_bcftools(sample_index, genome, run_info, reads_filename, genome_absolute_path, sample_platform, results_filename, soft_paths_dict, threads):
    if "/" in results_filename:
        print("Error: results_filename contains '/'", file=sys.stderr)
        sys.exit(1)
    tmp_root = get_tmp_root()
    tmp_bcftools_dir = os.path.join(tmp_root, sample_index, genome, run_info, "tmp_bcftools_results")
    if not os.path.exists(tmp_bcftools_dir):
        os.makedirs(tmp_bcftools_dir)
    if os.path.exists(results_filename):
        os.remove(results_filename)
    read_records = SeqIO.parse(reads_filename, "fasta")
    with cf.ThreadPoolExecutor(int(threads)) as tex:
        futures = [tex.submit(run_bcftools, read_record, genome_absolute_path, sample_platform, soft_paths_dict, tmp_bcftools_dir) for read_record in read_records]
        results = [future.result() for future in cf.as_completed(futures)]
    tmp_results_file = os.path.join(tmp_bcftools_dir, results_filename)
    if not os.path.exists(tmp_results_file):
        open(tmp_results_file, "wt").close()
    command_1 = "mv " + tmp_bcftools_dir + "/" + results_filename + " ."
    command_2 = "rm -rf " + os.path.join(tmp_root, sample_index)
    commands = command_1 + ";" + command_2
    hfbase.get_cli_output_lines(commands, side_effect = True)


def snv_or_indel(input_file, output_file):
    if os.path.exists(output_file):
        os.remove(output_file)
    output_columns = ["POS", "REF", "ALT", "TYPE", "ID_list", "counts"]
    if not os.path.exists(input_file) or os.path.getsize(input_file) == 0:
        df_reform = pl.DataFrame({column: [] for column in output_columns})
        df_reform.write_excel(output_file)
        return df_reform
    df = pl.read_csv(input_file, separator="\t", has_header=False, new_columns=["ID", "POS", "REF", "ALT"])
    df = df.with_columns(
        pl.when((pl.col("REF").str.len_chars() == 1) & (pl.col("ALT").str.len_chars() == 1))
        .then(pl.lit("SNV"))
        .otherwise(pl.lit("InDel"))
        .alias("TYPE")
    )
    df_reform = (
        df.group_by(["POS", "REF", "ALT", "TYPE"])
        .agg([
            pl.col("ID").alias("ID_list"),
            pl.count("ID").alias("counts")
        ])
        .with_columns(pl.col("ID_list").list.join(",").alias("ID_list"))
    )
    df_reform.write_excel(output_file)
    return df_reform
