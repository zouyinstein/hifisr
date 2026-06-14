# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

from Bio import SeqIO
import hifisr_functions.base as hfbase
import hifisr_functions.variants as hfvar
import numpy as np
from collections import OrderedDict
import sys
import os
import shutil
import time

# Function purity marker. "pure" means deterministic from explicit inputs with
# no file, shell, environment, logging, or input-mutation side effects.
FUNCTION_PURITY = {
    "replace_fasta_id": "impure",
    "get_subseq": "pure",
    "get_rc": "impure",
    "rotate_ref_to_non_repeat_region": "impure",
    "find_continous_zeros": "impure",
    "rotate_fasta": "impure",
    "mecat_cns": "impure",
    "flye_assemble": "impure",
    "simple_draft_asm": "impure",
    "flye_polish": "impure",
    "aln_to_ref": "impure",
    "get_flipped_fasta": "impure",
}


def replace_fasta_id(genome, input_fasta_path, output_fasta_path):
    records = list(SeqIO.parse(input_fasta_path, "fasta"))
    count = 1
    fout = open(output_fasta_path, "wt")
    for record in records:
        ID = genome + "_" + str(count) + " [" + record.description + "]"
        count += 1
        print(">" + ID, file=fout)
        print(record.seq, file=fout)
    if count >= 3:
        print("There are", count - 1, "contigs in the input reference.", file=sys.stderr)
    elif count == 2:
        print("There is", count - 1, "contig in the input reference.", file=sys.stderr)
    fout.close()


def get_subseq(ref, start, end, flank=0):
    subseq = ref[start-1:end]
    subseq_flank_start = start-1-flank
    subseq_flank_end = end+flank
    if subseq_flank_start < 0:
        subseq_flank_start = 0
    if subseq_flank_end > len(ref):
        subseq_flank_end = len(ref)
    subseq_flank = ref[subseq_flank_start:subseq_flank_end]
    return subseq, subseq_flank


def get_rc(input_fa_path, rc_fa_path):
    record = SeqIO.read(input_fa_path, "fasta")
    rc_record = record.reverse_complement(id=True, name=True, description=True)
    with open(rc_fa_path, "wt") as fout:
        print(rc_record.format("fasta"), file=fout)


def _write_lines(path, lines):
    with open(path, "wt") as fout:
        for line in lines:
            print(line, file=fout)


def _sorted_interval(start, end):
    left = min(int(start), int(end))
    right = max(int(start), int(end))
    return left, right


def _is_full_self_hit(fields, ref_len):
    q_start, q_end = _sorted_interval(fields[6], fields[7])
    s_start, s_end = _sorted_interval(fields[8], fields[9])
    aln_len = int(fields[3])
    pident = float(fields[2])
    return (
        q_start == 1
        and q_end == ref_len
        and s_start == 1
        and s_end == ref_len
        and aln_len >= ref_len
        and pident >= 99.999
    )


def _build_repeat_mask_from_blastn(blastn_lines, ref_len):
    diff = np.zeros(ref_len + 1, dtype=int)
    used_alignment_count = 0
    skipped_full_self_count = 0
    used_rows = []

    for line in blastn_lines:
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        if _is_full_self_hit(fields, ref_len):
            skipped_full_self_count += 1
            continue
        q_start, q_end = _sorted_interval(fields[6], fields[7])
        start = max(q_start - 1, 0)
        end = min(q_end, ref_len)
        if start >= end:
            continue
        diff[start] += 1
        diff[end] -= 1
        used_alignment_count += 1
        used_rows.append(line)

    repeat_mask = np.cumsum(diff[:-1])
    stats = {
        "used_alignment_count": used_alignment_count,
        "skipped_full_self_count": skipped_full_self_count,
    }
    return repeat_mask, stats, used_rows


def _intervals_from_mask(mask, want_repeat):
    intervals = []
    start = None
    for index, value in enumerate(mask):
        is_target = value > 0 if want_repeat else value == 0
        if is_target and start is None:
            start = index
        elif not is_target and start is not None:
            intervals.append((start, index - 1, index - start))
            start = None
    if start is not None:
        intervals.append((start, len(mask) - 1, len(mask) - start))
    return intervals


def _circular_non_repeat_blocks(non_repeat_intervals, ref_len):
    if not non_repeat_intervals:
        return []
    blocks = list(non_repeat_intervals)
    if len(blocks) > 1 and blocks[0][0] == 0 and blocks[-1][1] == ref_len - 1:
        first = blocks.pop(0)
        last = blocks.pop(-1)
        blocks.append((last[0], first[1], last[2] + first[2]))
    blocks.sort(key=lambda item: item[2], reverse=True)
    return blocks


def _midpoint_on_circle(block, ref_len):
    start, end, _length = block
    if start <= end:
        return (start + end) // 2
    return ((start + end + ref_len) // 2) % ref_len


def _write_bed(path, record_id, intervals, label):
    with open(path, "wt") as fout:
        for start, end, length in intervals:
            print(record_id, start, end + 1, label, length, sep="\t", file=fout)


def _write_reference_rotation_summary(path, summary):
    keys = list(summary)
    with open(path, "wt") as fout:
        print("\t".join(keys), file=fout)
        print("\t".join(str(summary[key]) for key in keys), file=fout)


def rotate_ref_to_non_repeat_region(genome, genome_fasta_path, soft_paths_dict, rotation=False):
    ref_records = list(SeqIO.parse(genome_fasta_path, "fasta"))
    if len(ref_records) != 1:
        print("Error: fasta has more than one record.")
        return
    ref_len = len(ref_records[0].seq)
    commands = soft_paths_dict.get("blastn") + " -query " + genome_fasta_path + " -subject " + genome_fasta_path + " -outfmt 6"
    blastn_start = time.monotonic()
    blastn_lines = hfbase.get_cli_output_lines(commands, side_effect = False)
    blastn_seconds = time.monotonic() - blastn_start
    _write_lines("blastn_self.tsv", blastn_lines)
    repeat_pos_array, stats, used_rows = _build_repeat_mask_from_blastn(blastn_lines, ref_len)
    _write_lines("blastn_self_used_for_repeat_mask.tsv", used_rows)
    repeat_intervals = _intervals_from_mask(repeat_pos_array, want_repeat=True)
    non_repeat_intervals = _intervals_from_mask(repeat_pos_array, want_repeat=False)
    circular_blocks = _circular_non_repeat_blocks(non_repeat_intervals, ref_len)
    _write_bed("repeat_regions.bed", ref_records[0].id, repeat_intervals, "repeat")
    _write_bed("non_repeat_regions.linear.bed", ref_records[0].id, non_repeat_intervals, "non_repeat")
    _write_bed("non_repeat_regions.circular_ranked.bed", ref_records[0].id, circular_blocks, "non_repeat_circular")

    if not circular_blocks:
        print("Error: no non-repeat region found.", file=sys.stderr)
        return

    non_start = circular_blocks[0][0]
    non_end = circular_blocks[0][1]
    length = circular_blocks[0][2]
    rot_step = _midpoint_on_circle(circular_blocks[0], ref_len)
    summary = {
        "genome": genome,
        "method": "blastn",
        "blastn_options": "-outfmt 6",
        "min_non_repeat_len": 5000,
        "reference_length": ref_len,
        "blastn_seconds": f"{blastn_seconds:.6f}",
        "blastn_alignment_count": len(blastn_lines),
        "used_alignment_count": stats["used_alignment_count"],
        "skipped_full_self_count": stats["skipped_full_self_count"],
        "repeat_covered_bp": int(np.count_nonzero(repeat_pos_array)),
        "non_repeat_bp": int(ref_len - np.count_nonzero(repeat_pos_array)),
        "repeat_region_count": len(repeat_intervals),
        "non_repeat_region_count_linear": len(non_repeat_intervals),
        "non_repeat_region_count_circular": len(circular_blocks),
        "longest_non_repeat_start_0based": non_start,
        "longest_non_repeat_end_0based": non_end,
        "longest_non_repeat_length": length,
        "rot_step": rot_step,
        "rotated_fasta": genome + "_rotated_" + str(rot_step) + ".fasta" if length > 5000 else "",
    }
    _write_reference_rotation_summary("reference_rotation_summary.tsv", summary)

    if rotation:
        if length > 5000:
            print("The largest non-repeat region is " + str(length) + " bp long.", file=sys.stderr)
            rotate_fasta(genome_fasta_path, genome + "_rotated_" + str(rot_step) + ".fasta", rot_step)
        else:
            print("The largest non-repeat region is less than 5 kb.", file=sys.stderr)
        return rot_step
    else:
        return 0


def find_continous_zeros(info_list, repeat_pos_array):
    for i in range(0, len(repeat_pos_array)):
        if repeat_pos_array[i] == 0:
            piece_start = i
            piece_end = i
            repeat_pos_array[i] = 1
            break
    for i in range(piece_start+1, len(repeat_pos_array)):
        if repeat_pos_array[i] == 0:
            piece_end = i
            repeat_pos_array[i] = 1
            if i == len(repeat_pos_array) - 1:
                length = piece_end - piece_start + 1
                info_list.append((piece_start, piece_end, length))
                break
        else:
            length = piece_end - piece_start + 1
            info_list.append((piece_start, piece_end, length))
            break
    return info_list, repeat_pos_array


def rotate_fasta(genome_fasta_path, rotated_fasta_path, step):
    ref_records = list(SeqIO.parse(genome_fasta_path, "fasta"))
    if len(ref_records) != 1:
        print("Error: input fasta has more than one record.")
        return
    id = ref_records[0].id
    sequence = ref_records[0].seq
    if step > 0 and step < len(sequence):
        subseq_1 = sequence[step:]
        subseq_2 = sequence[0:step]
        rot_seq = subseq_1 + subseq_2
    elif step == 0:
        rot_seq = sequence
    elif step < 0:
        subseq_1 = sequence[:(len(sequence)+step)]
        subseq_2 = sequence[(len(sequence)+step):]
        rot_seq = subseq_2 + subseq_1
    else:
        print("Error: step must be greater than 0 and less than the length of the sequence.")
        return None
    with open(rotated_fasta_path, "wt") as fout:
        print(">" + id + " [rotation=" + str(step) + "]", file=fout)
        print(rot_seq, file=fout)


def mecat_cns(genome, genome_size, reads, soft_paths_dict, threads):
    command_0 = "rm -rf mecat_" + genome + "_" + str(genome_size) + " " + genome + "_" + str(genome_size) + "_config.txt"
    command_1 = soft_paths_dict.get("mecat") + " config " + genome + "_" + str(genome_size) + "_config.txt"
    command_2 = "sed s/PROJECT=/PROJECT=mecat_" + genome + "_" + str(genome_size) + "/ " + genome + "_" + str(genome_size) + "_config.txt > new_config && mv new_config " + genome + "_" + str(genome_size) + "_config.txt"
    command_3 = "sed s/RAWREADS=/RAWREADS=" + reads + "/ " + genome + "_" + str(genome_size) + "_config.txt > new_config && mv new_config " + genome + "_" + str(genome_size) + "_config.txt"
    command_4 = "sed s/GENOME_SIZE=/GENOME_SIZE=" + str(genome_size) + "000/ " + genome + "_" + str(genome_size) + "_config.txt > new_config && mv new_config " + genome + "_" + str(genome_size) + "_config.txt"
    command_5 = "sed s/THREADS=/THREADS=" + threads + "/ " + genome + "_" + str(genome_size) + "_config.txt > new_config && mv new_config " + genome + "_" + str(genome_size) + "_config.txt"
    command_6 = soft_paths_dict.get("mecat") + " correct " + genome + "_" + str(genome_size) + "_config.txt"
    command_7 = "cp mecat_" + genome + "_" + str(genome_size) + "/1-consensus/cns_final.fasta mecat_" + genome + "_" + str(genome_size) + ".fasta"
    command_8 = "rm -rf mecat_" + genome + "_" + str(genome_size) + " mecat_" + genome + "_" + str(genome_size) + "_config.txt"
    commands = command_0 + " && " + command_1 + " && " + command_2 + " && " + command_3 + " && " + command_4 + " && " + command_5 + " && " + command_6 + " && " + command_7 + " && " + command_8
    hfbase.run_checked(commands)
    return


def flye_assemble(prefix, genome, genome_size, reads, soft_paths_dict, sample_platform, threads, correction=True):
    if correction:
        platform_dict = {
            "HiFi": "--pacbio-hifi",
            "CLR": "--pacbio-corr",
            "ONT": "--nano-corr",
            "ultra-long": "--nano-corr",
        }
    else:
        platform_dict = {
            "HiFi": "--pacbio-hifi",
            "CLR": "--pacbio-raw",
            "ONT": "--nano-raw",
            "ultra-long": "--nano-raw",
        }
    flye = soft_paths_dict.get("flye", "flye")
    # command_1 = flye + " --meta " + platform_dict.get(sample_platform) + " " + reads + " --extra-params output_gfa_before_rr=1 --genome-size " + str(genome_size) + "K -t " + threads + " -o flye_" + prefix + "_" + genome + "_" + str(genome_size) + "K"
    command_1 = flye + " " + platform_dict.get(sample_platform) + " " + reads + " --extra-params output_gfa_before_rr=1 --genome-size " + str(genome_size) + "K -t " + threads + " -o flye_" + prefix + "_" + genome + "_" + str(genome_size) + "K"
    command_2 = "cp flye_" + prefix + "_" + genome + "_" + str(genome_size) + "K/assembly_graph.gfa " + prefix + "_" + genome + "_" + str(genome_size) + "K_after_rr.gfa"
    command_3 = "cp flye_" + prefix + "_" + genome + "_" + str(genome_size) + "K/20-repeat/graph_before_rr.gfa " + prefix + "_" + genome + "_" + str(genome_size) + "K_before_rr.gfa"
    command_4 = "rm -rf flye_" + prefix + "_" + genome + "_" + str(genome_size) + "K"
    commands = command_1 + " && " + command_2 + " && " + command_3 + " && " + command_4
    hfbase.run_checked(commands)
    return


def simple_draft_asm(genome, genome_size, preset, reads, soft_paths_dict, threads):
    output_prefix = "simple_draft_asm_" + genome + "_" + str(genome_size) + "K_" + preset
    output_dir = output_prefix + "_work"
    output_gfa = output_prefix + ".gfa"
    simple_draft_asm_bin = soft_paths_dict.get("simple_draft_asm", "simple_draft_asm")
    hfbase.run_checked("rm -rf " + output_dir + " " + output_gfa)
    command_1 = simple_draft_asm_bin + " -p " + preset + " -i " + reads + " -o " + output_dir + " -t " + threads
    hfbase.run_checked(command_1)
    graph_gfa = os.path.join(output_dir, "graph.gfa")
    if not os.path.exists(graph_gfa):
        subset_graphs = []
        for name in os.listdir(output_dir):
            if not name.startswith("read_subset_"):
                continue
            subset_graph = os.path.join(output_dir, name, "graph.gfa")
            if os.path.exists(subset_graph):
                try:
                    subset_size = int(name.rsplit("_", 1)[1])
                except ValueError:
                    subset_size = 0
                subset_graphs.append((subset_size, subset_graph))
        if subset_graphs:
            graph_gfa = sorted(subset_graphs)[-1][1]
    if not os.path.exists(graph_gfa):
        raise FileNotFoundError("simple_draft_asm did not produce graph.gfa for preset " + preset)
    shutil.copyfile(graph_gfa, output_gfa)
    hfbase.run_checked("rm -rf " + output_dir)
    return output_gfa


def flye_polish(genome, before_fasta_path, after_fasta_prefix, reads, soft_paths_dict, sample_platform, threads, correction=True):
    if correction:
        platform_dict = {
            "HiFi": "--pacbio-hifi",
            "CLR": "--pacbio-corr",
            "ONT": "--nano-corr",
            "ultra-long": "--nano-corr",
        }
    else:
        platform_dict = {
            "HiFi": "--pacbio-hifi",
            "CLR": "--pacbio-raw",
            "ONT": "--nano-raw",
            "ultra-long": "--nano-raw",
        }
    after_fasta_path = after_fasta_prefix + ".fasta"
    flye = soft_paths_dict.get("flye", "flye")
    command_1 = flye + " --meta " + platform_dict.get(sample_platform) + " " + reads + " --polish-target " + before_fasta_path + " -i 2 -t " + threads + " -o " + genome + "_flye_polish"
    command_2 = "cp " + genome + "_flye_polish/polished_2.fasta " + after_fasta_path + " && rm -rf " + genome + "_flye_polish"
    command_3 = soft_paths_dict.get("minimap2") + " -t 1 -ax map-hifi " + before_fasta_path + " " + after_fasta_path
    command_4 = soft_paths_dict.get("samtools") + " view -Sb -F 4 -@ 1 -"
    command_5 = soft_paths_dict.get("samtools") + " sort -@ 1 -o -"
    command_6 = soft_paths_dict.get("bcftools") + " mpileup --indels-2.0 -m 1 -Ou -f " + before_fasta_path + " -"
    command_7 = soft_paths_dict.get("bcftools") + " call -mv -P 0.99 -Ov | grep -v '^#' | cut -f2,4,5 > " + after_fasta_prefix + "_variants.txt"
    commands = command_1 + " && " + command_2 + " && " + command_3 + " | " + command_4 + " | " + command_5 + " | " + command_6 + " | " + command_7
    ret = hfbase.get_cli_output_lines(commands, side_effect = True)
    return


def aln_to_ref(genome, genome_absolute_path, polished_fasta_absolute_path, final_fasta_name, results_filename, soft_paths_dict):
    tmp_blastn_dir = "tmp_blastn_results"
    if os.path.exists(tmp_blastn_dir):
        hfbase.get_cli_output_lines("rm -rf " + tmp_blastn_dir)
    os.makedirs(tmp_blastn_dir)
    if "/" in results_filename:
        print("Error: results_filename contains '/'", file=sys.stderr)
        sys.exit(1)
    if os.path.exists(results_filename):
        os.remove(results_filename)    
    replace_fasta_id(genome, genome_absolute_path, "ref_" + genome + ".fa")
    replace_fasta_id(genome, polished_fasta_absolute_path, "before_" + genome + ".fa")
    ref_record = SeqIO.read("ref_" + genome + ".fa", "fasta")
    hfvar.run_blastn_sorter_single(ref_record, "before_" + genome + ".fa", results_filename, soft_paths_dict, tmp_blastn_dir)
    command_1 = "mv " + tmp_blastn_dir + "/" + results_filename + " ."
    command_2 = "rm -rf " + tmp_blastn_dir
    commands = command_1 + ";" + command_2
    hfbase.get_cli_output_lines(commands, side_effect = True)

    adj_count = 0
    rc_or_not = 0
    rot_or_not = 0
    # rc or not
    blastn_results = hfbase.get_file_lines(results_filename)
    fields = blastn_results[0].split("\t")
    aln_1_info = fields[5].split(";")
    tmp_od = OrderedDict()
    for tmp_aln_info in aln_1_info:
        tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
        tmp_od[tmp_aln_info_key] = tmp_aln_info_value
    if tmp_od.get("strand") == "-":
        adj_count += 1
        rc_or_not = 1
        get_rc("before_" + genome + ".fa", "before_" + genome + "_rc.fa")
        replace_fasta_id(genome + "_adj_" + str(adj_count), "before_" + genome + "_rc.fa", "adj_" + str(adj_count) + ".fa")
        hfbase.get_cli_output_lines("cp adj_" + str(adj_count) + ".fa adjust.fa", side_effect = True)
        if os.path.exists(results_filename):
            os.remove(results_filename)
        if os.path.exists(tmp_blastn_dir):
            hfbase.get_cli_output_lines("rm -rf " + tmp_blastn_dir)
        os.makedirs(tmp_blastn_dir)
        hfvar.run_blastn_sorter_single(ref_record, "adjust.fa", results_filename, soft_paths_dict, tmp_blastn_dir)
        command_1 = "mv " + tmp_blastn_dir + "/" + results_filename + " ."
        command_2 = "rm -rf " + tmp_blastn_dir
        commands = command_1 + ";" + command_2
        hfbase.get_cli_output_lines(commands, side_effect = True)

    # rotate or not
    blastn_results = hfbase.get_file_lines(results_filename)
    fields = blastn_results[0].split("\t")
    aln_1_info = fields[5].split(";")
    tmp_od = OrderedDict()
    for tmp_aln_info in aln_1_info:
        tmp_aln_info_key, tmp_aln_info_value = tmp_aln_info.split("=")
        tmp_od[tmp_aln_info_key] = tmp_aln_info_value
    if tmp_od.get("strand") == "-":
        print("Error: the first alignment is on the reverse strand.")
        sys.exit(1)
    qs = int(tmp_od.get("qs"))
    ss = int(tmp_od.get("ss"))
    rot_step = int(ss - qs)
    if rot_step != 0:
        adj_count += 1
        rot_or_not = 1
        if adj_count == 1:
            rotate_fasta("before_" + genome + ".fa", "adj_" + str(adj_count) + ".fa", rot_step)
        else:
            rotate_fasta("adjust.fa", "adj_" + str(adj_count) + ".fa", rot_step)
        replace_fasta_id(genome + "_adj_" + str(adj_count), "adj_" + str(adj_count) + ".fa", "adj_" + str(adj_count) + ".fa")
        hfbase.get_cli_output_lines("cp adj_" + str(adj_count) + ".fa adjust.fa", side_effect = True)
        if os.path.exists(results_filename):
            os.remove(results_filename)
        if os.path.exists(tmp_blastn_dir):
            hfbase.get_cli_output_lines("rm -rf " + tmp_blastn_dir)
        os.makedirs(tmp_blastn_dir)
        hfvar.run_blastn_sorter_single(ref_record, "adjust.fa", results_filename, soft_paths_dict, tmp_blastn_dir)
        command_1 = "mv " + tmp_blastn_dir + "/" + results_filename + " ."
        command_2 = "rm -rf " + tmp_blastn_dir
        commands = command_1 + ";" + command_2
        hfbase.get_cli_output_lines(commands, side_effect = True)

    if adj_count == 0:
        hfbase.get_cli_output_lines("mv before_" + genome + ".fa " + genome + "_aligned.fasta", side_effect = True)
        hfbase.get_cli_output_lines("rm ref_" + genome + ".fa", side_effect = True)
    else:
        hfbase.get_cli_output_lines("rm before_" + genome + ".fa ref_" + genome + ".fa adj_*.fa", side_effect = True)
        hfbase.get_cli_output_lines("mv adjust.fa " + genome + "_aligned.fasta", side_effect = True)
    command_1 = "mv " + genome + "_aligned.fasta " + final_fasta_name
    hfbase.get_cli_output_lines(command_1, side_effect = True)
    return adj_count, rc_or_not, rot_or_not


def get_flipped_fasta(genome_fasta_path, flipped_fasta_path, flip_start, flip_end):
    flip_start = int(flip_start)
    flip_end = int(flip_end)
    record = SeqIO.read(genome_fasta_path, "fasta")
    id = record.id
    sequence = record.seq
    if flip_start >= 1 and flip_end <= len(sequence):
        subseq_1 = sequence[0:flip_start]
        subseq_2 = sequence[flip_start:(flip_end-1)]
        subseq_3 = sequence[(flip_end-1):]
        flipped_seq = subseq_1 + subseq_2.reverse_complement() + subseq_3
        with open(flipped_fasta_path, "wt") as fout:
            print(">" + id + " [flip=" + str(flip_start) + "-" + str(flip_end) + "]", file=fout)
            print(flipped_seq, file=fout)
        print("The fasta file has been flipped via coordinates " + str(flip_start) + " and " + str(flip_end) + ".", file=sys.stderr)
    return
