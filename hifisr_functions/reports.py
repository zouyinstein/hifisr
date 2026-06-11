# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

import hifisr_functions.base as hfbase
import numpy as np
import pandas as pd
from Bio import SeqIO
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "hifisr_matplotlib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from collections import OrderedDict
import sys

# Function purity marker. "pure" means deterministic from explicit inputs with
# no file, shell, environment, logging, or input-mutation side effects.
FUNCTION_PURITY = {
    "get_fastq_stats": "impure",
    "plot_length_qual": "impure",
    "plot_bubble_type_2_rep_raw": "impure",
    "plot_coverage": "impure",
    "get_gfa_blastn_png": "impure",
    "convert_blastn_alignments_to_table": "impure",
}


def get_fastq_stats(prefix, sample_fastq_path, soft_paths_dict, threads):
    commands = soft_paths_dict.get("seqkit") + " stat -a -T " + sample_fastq_path + " -j " + threads
    output_lines_1 = hfbase.get_cli_output_lines(commands)
    total_read_number = output_lines_1[1].split("\t")[3]
    total_bases = output_lines_1[1].split("\t")[4]
    id_length_qual_file = prefix + "_id_length_qual.txt"
    commands = soft_paths_dict.get("seqkit") + " fx2tab -j " + threads + " -qlni " + sample_fastq_path + " > " + id_length_qual_file
    hfbase.get_cli_output_lines(commands) # ID, length, qual 
    return id_length_qual_file, total_read_number, total_bases


def plot_length_qual(prefix, sample_platform, id_length_qual_file, total_read_number, total_bases):
    id_length_qual_lines = hfbase.get_file_lines(id_length_qual_file)
    length_list = [ int(line.split("\t")[1]) for line in id_length_qual_lines ]
    qual_list = [ float(line.split("\t")[2]) for line in id_length_qual_lines ]
    length_bin_dict = {
        "HiFi": (50000, 500),
        "CLR": (50000, 500),
        "ONT": (50000, 500),
        "ultra-long": (100000, 500),
        "Short": (500, 10),
    }
    fig = plt.figure(figsize=(10, 8))
    gs = gridspec.GridSpec(13, 2)
    ax1 = plt.subplot(gs[0:5, :])
    ax2 = plt.subplot(gs[7:12, :])
    axs = [ax1, ax2]
    # plot length distribution
    length_bins = np.arange(0, length_bin_dict.get(sample_platform)[0], length_bin_dict.get(sample_platform)[1])
    axs[0].hist(length_list, bins=length_bins, color="blue", alpha=0.5, edgecolor="black", linewidth=0.5)
    axs[0].axvline(sum(length_list)/len(length_list), color="red", linestyle="--")
    axs[0].set_xlim([0, length_bin_dict.get(sample_platform)[0]])
    axs[0].grid(True, alpha=0.5)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    axs[0].set_xlabel("Read Length")
    axs[0].set_ylabel("Counts")
    axs[0].set_title("total_read_number = " + total_read_number + "; total_bases = " + total_bases)
    axs[1].hist(qual_list, bins=np.arange(0, 100, 1), color="red", alpha=0.5, edgecolor="black", linewidth=0.5)
    axs[1].axvline(sum(qual_list)/len(qual_list), color="red", linestyle="--")
    axs[1].set_xlim([0, 100])
    axs[1].grid(True, alpha=0.5)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    axs[1].set_xlabel("Read Quality")
    axs[1].set_ylabel("Counts")
    axs[1].set_title("total_read_number = " + total_read_number + "; total_bases = " + total_bases)
    plt.savefig(prefix + "_length_qual_distribution.pdf")
    # clear the figure
    plt.clf()
    # create 2-D histogram for length and quality
    fig = plt.figure(figsize=(10, 8))
    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[0, 0])
    white_to_red = mcolors.LinearSegmentedColormap.from_list("white_red", ["white", "red"])
    hb = plt.hexbin(length_list, qual_list, gridsize=30, cmap=white_to_red, edgecolors="lightgray")
    cb = plt.colorbar(hb, label="Count in bin")
    ax.set_xlabel("Read Length")
    ax.set_ylabel("Read Quality")
    ax.set_title("total_read_number = " + total_read_number + "; total_bases = " + total_bases)
    plt.savefig(prefix + "_length_qual_2d_distribution.pdf")
    plt.close()
    return


def plot_bubble_type_2_rep_raw(table_file, IDs_dir, ref_fasta):
    if os.path.exists(table_file):
        df = pd.read_excel(table_file, sheet_name='Sheet1') 
    else:
        return
    FL_ids_files = [IDs_dir + "/" + file for file in os.listdir(IDs_dir) if file.endswith("_FL_ids.txt")]
    FL_count = 0
    for file in FL_ids_files:
        FL_count += len(hfbase.get_file_lines(file))
    ref_len = len(SeqIO.read(ref_fasta, "fasta").seq) # mito_rotated_flye_polish_1.fasta
    df['subgroup_count_norm'] = df['subgroup_count'] / FL_count * 10000  # normalize to 10000 FL reads
    fig = plt.figure(figsize=(10, 10)) # set figure size with
    df["color"] = df["mid_olp_1"].apply(lambda x: "#FF00FF" if x >= 1000 else "#00FF00" if x >= 300 else "#0016FF" if x >= 200 else "#E8720C" if x >= 100 else "#E6E600" if x >= 50 else "#B0B0B0") # set colors for the bubble plot by the mid_olp_1
    df["se1"] = df["(se1, ss2)"].apply(lambda x: int(x.split(",")[0].lstrip("("))) # (se1, ss2)
    df["ss2"] = df["(se1, ss2)"].apply(lambda x: int(x.split(",")[1].rstrip(")"))) # (se1, ss2)
    size_ratio = 10
    alpha_value = 0.5 # df["alpha"] = 0.5
    if len(df) > 0:
        plt.scatter(df['se1'], df['ss2'], s=df['subgroup_count_norm']*size_ratio, c=df['color'], alpha=alpha_value, linewidths=0)
        plt.grid(True, alpha=0.5)
        plt.xlim(1, ref_len)
        plt.ylim(1, ref_len)
        plt.savefig('bubble_type_2_rep_raw.pdf')
    plt.close()
    return


def plot_coverage(cov_file_1, cov_file_2, cov_file_3, start, end, fig_length=12, fig_height=3):
    if os.path.exists(cov_file_1):
        cov_1 = np.loadtxt(cov_file_1, dtype=int, usecols=1)
    else:
        print("No coverage file found for " + cov_file_1, file=sys.stderr)
    if os.path.exists(cov_file_2) and os.path.getsize(cov_file_2) > 0: # and size not zero
        cov_2 = np.loadtxt(cov_file_2, dtype=int, usecols=1)
    else:
        cov_2 = np.zeros(cov_1.shape)
    if os.path.exists(cov_file_3) and os.path.getsize(cov_file_3) > 0: # and size not zero
        cov_3 = np.loadtxt(cov_file_3, dtype=int, usecols=1)
    else:
        cov_3 = np.zeros(cov_1.shape)
    cov_combine = cov_1 + cov_2
    fig = plt.figure(figsize=(fig_length, fig_height), dpi=600)
    plt.plot(cov_combine[(start-1):end], color="#EAB13E", label="FL")
    plt.plot(cov_1[(start-1):end], color="#D1D1D1", label="partial")
    plt.fill_between(np.arange(start-1, end), cov_combine[(start-1):end], 1, color="#EAB13E", alpha=1) # run long time
    plt.fill_between(np.arange(start-1, end), cov_1[(start-1):end], 1, color="#D1D1D1", alpha=1)
    plt.plot(cov_3[(start-1):end], color="#5CAB38", label="variant", linewidth=1)
    plt.grid(True, alpha=0.5)
    ax = plt.gca()
    ax.set_xlim([start, end+1])
    max_y = int(np.max(cov_combine[(start-1):end]))
    ax.set_ylim([0, max_y+100])
    plt.savefig('coverage_' + str(start) + "_" + str(end) + '.pdf')
    plt.savefig('coverage_' + str(start) + "_" + str(end) + '.png')
    plt.close()
    return


def get_gfa_blastn_png(genome_absolute_path, soft_paths_dict):
    command_3 = "ls -1 *.gfa | while read i; do " + soft_paths_dict.get("bandage") + " image $i ${i%.gfa}.png --edgelen 25 --singlearr --depwidth 0.8 --colour blastsolid --query " + genome_absolute_path + "; done"
    hfbase.run_checked(command_3)
    return


def convert_blastn_alignments_to_table(blastn_alignments_file, output_file):
    # refine into Polars in the future
    for alignment_line in hfbase.get_file_lines(blastn_alignments_file):
        df_blastn_alignments = pd.DataFrame(columns=["aln_index", "aln_len", "aln_olp_len", "aln_idt", "aln_strand", "aln_qs", "aln_qe", "aln_ss", "aln_se", "aln_cn"])
        _, _, _, aln_type, percent_total, *align_info_list_of_dict = alignment_line.split("\t")
        for i in range(len(align_info_list_of_dict)):
            align_info_list_of_dict[i] = OrderedDict([ (x.split("=")[0], x.split("=")[1]) for x in align_info_list_of_dict[i].split(";") ])
            # aln=1;len=49386;olp=26264;idt=99.872;strand=+;qs=1;qe=49386;ss=1;se=49397;cn=1;c1=100.0,1,49397
            aln_index = align_info_list_of_dict[i]["aln"] # a string
            aln_len = align_info_list_of_dict[i]["len"]
            aln_olp_len = align_info_list_of_dict[i]["olp"]
            aln_idt = align_info_list_of_dict[i]["idt"]
            aln_strand = align_info_list_of_dict[i]["strand"]
            aln_qs = align_info_list_of_dict[i]["qs"]
            aln_qe = align_info_list_of_dict[i]["qe"]
            aln_ss = align_info_list_of_dict[i]["ss"]
            aln_se = align_info_list_of_dict[i]["se"]
            aln_cn = align_info_list_of_dict[i]["cn"]
            df_blastn_alignments.loc[i] = [aln_index, aln_len, aln_olp_len, aln_idt, aln_strand, aln_qs, aln_qe, aln_ss, aln_se, aln_cn]
    df_blastn_alignments.to_excel(output_file, index=False) # "blastn_alignments.xlsx"
    return
