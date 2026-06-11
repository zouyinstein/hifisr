# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

import hifisr_functions.base as hfbase
import hifisr_functions.variants as hfvar
import hifisr_functions.references as hfref
import concurrent.futures as cf
from collections import OrderedDict
import pandas as pd
import polars as pl
import numpy as np
from Bio import SeqIO
import math
import os
import sys

# Function purity marker. "pure" means deterministic from explicit inputs with
# no file, shell, environment, logging, or input-mutation side effects.
FUNCTION_PURITY = {
    "get_variant_types": "impure",
    "combine_variant_anno": "impure",
    "add_depth_and_frq": "impure",
    "run_transfer_blastn": "impure",
    "merge_fragments_files": "impure",
    "remove_inner_fragments": "impure",
    "merge_numt_nupt": "impure",
    "get_raw_transfer_groups": "impure",
}


def get_variant_types(ref_fasta_file, input_excel, output_excel):
    ref_record = SeqIO.read(ref_fasta_file, "fasta")
    df_final_variants = pd.read_excel(input_excel)
    if len(df_final_variants) == 0:
        for column in ["annotation", "trinucleotide_context", "SNV_group", "InDel_group", "summary_anno"]:
            df_final_variants[column] = pd.Series(dtype="object")
        df_final_variants.to_excel(output_excel, index=False)
        return
    for i in range(len(df_final_variants)):
        pos = int(df_final_variants.loc[i, "POS"])
        ref = df_final_variants.loc[i, "REF"]
        alt = df_final_variants.loc[i, "ALT"]
        variant_type = df_final_variants.loc[i, "TYPE"]
        df_final_variants.loc[i, "annotation"] = "No"
        if variant_type == "SNV":
            if pos == 1:
                df_final_variants.loc[i, "trinucleotide_context"] = str(ref_record.seq[-1] + ref + ref_record.seq[pos])
            elif pos == len(ref_record.seq):
                df_final_variants.loc[i, "trinucleotide_context"] = str(ref_record.seq[pos-2] + ref + ref_record.seq[0])
            else:
                df_final_variants.loc[i, "trinucleotide_context"] = str(ref_record.seq[pos-2] + ref + ref_record.seq[pos])
            df_final_variants.loc[i, "SNV_group"] = str(ref) + ">" + str(alt)
            df_final_variants.loc[i, "InDel_group"] = "-"
            df_final_variants.loc[i, "summary_anno"] = "-"
            df_final_variants.loc[i, "annotation"] = "Yes"
        if variant_type == "InDel":
            df_final_variants.loc[i, "trinucleotide_context"] = "-"
            df_final_variants.loc[i, "SNV_group"] = "-"
            df_final_variants.loc[i, "InDel_group"] = "-"
            df_final_variants.loc[i, "summary_anno"] = "-" # modify later
            # for homopolymer: unit_size = 1, len(ref) >1, ref_base = ref[-1], ref_size = len(ref) - 1, alt_size = len(alt) -1
            if len(ref) > 1 and len(alt) > 1:
                ref_str = ref[1:]
                alt_str = alt[1:]
                ref_unit_base = ref[1]
                alt_unit_base = alt[1]
                if ref_unit_base == alt_unit_base:
                    unit_base = ref_unit_base
                    # count of unit_base in ref and alt
                    ref_count = ref_str.count(unit_base)
                    alt_count = alt_str.count(unit_base)
                    if ref_count == len(ref_str) and alt_count == len(alt_str):
                        df_final_variants.loc[i, "InDel_group"] = "homopolymer"
                        df_final_variants.loc[i, "summary_anno"] = "poly-" + str(unit_base) + ";ref_size=" + str(len(ref_str)) + ";alt_size=" + str(len(alt_str))
                        df_final_variants.loc[i, "annotation"] = "Yes"
            # for homodimer
            if len(ref) > 2 and len(alt) > 2:
                ref_str = ref[1:]
                alt_str = alt[1:]
                ref_unit_base = ref[1:3]
                alt_unit_base = alt[1:3]
                if ref_unit_base == alt_unit_base and ref_unit_base[0] != ref_unit_base[1]: # avoid the same base
                    unit_base = ref_unit_base
                    # count of unit_base in ref and alt
                    ref_count = len(ref_str)//2
                    alt_count = len(alt_str)//2
                    if ref_str == unit_base * ref_count and alt_str == unit_base * alt_count:
                        df_final_variants.loc[i, "InDel_group"] = "homodimer"
                        df_final_variants.loc[i, "summary_anno"] = "dimer-" + str(unit_base) + ";ref_size=" + str(ref_count) + ";alt_size=" + str(alt_count)
                        df_final_variants.loc[i, "annotation"] = "Yes"
            
            # for other cases
            if df_final_variants.loc[i, "annotation"] == "No":
                # for tandem: unit_size > 1, ref_copy, alt_copy — slippage is here, add later
                # gain/loss of tandem copies for cases when one of them is a single copy unit
                if len(ref) > 1 and len(alt) > 1:
                    ref_str = ref[1:]
                    alt_str = alt[1:]
                    # guess the unit based on gcd
                    current_gcd = math.gcd(len(ref_str), len(alt_str)) # the most frequent count is 1
                    ref_count = len(ref_str) // current_gcd
                    alt_count = len(alt_str) // current_gcd
                    ref_unit = ref_str[:current_gcd]
                    if ref[0] == alt[0] and ref_str == ref_unit * ref_count and alt_str == ref_unit * alt_count:
                        df_final_variants.loc[i, "InDel_group"] = "tandem"
                        df_final_variants.loc[i, "summary_anno"] = "tandem-" + str(ref_unit) + ";ref_size=" + str(ref_count) + ";alt_size=" + str(alt_count)
                        df_final_variants.loc[i, "annotation"] = "Yes"
            if df_final_variants.loc[i, "annotation"] == "No":
                # also tamdem, but the first few bases are different, not 1 as above
                leading_base_count = 2
                while ref[:leading_base_count] == alt[:leading_base_count]:
                    if len(ref) > leading_base_count and len(alt) > leading_base_count:
                        ref_str = ref[leading_base_count:]
                        alt_str = alt[leading_base_count:]
                        current_gcd = math.gcd(len(ref_str), len(alt_str))
                        ref_count = len(ref_str) // current_gcd
                        alt_count = len(alt_str) // current_gcd
                        ref_unit = ref_str[:current_gcd]
                        if ref_str == ref_unit * ref_count and alt_str == ref_unit * alt_count:
                            df_final_variants.loc[i, "InDel_group"] = "tandem"
                            df_final_variants.loc[i, "summary_anno"] = "tandem-" + str(ref_unit) + ";ref_size=" + str(ref_count) + ";alt_size=" + str(alt_count)
                            df_final_variants.loc[i, "annotation"] = "Yes"
                            break
                    leading_base_count += 1
                    if leading_base_count == len(ref) or leading_base_count == len(alt):
                        break
            if df_final_variants.loc[i, "annotation"] == "No":
                leading_base_count = 1
                while ref[:leading_base_count] == alt[:leading_base_count]:
                    leading_base_count += 1
                leading_base_count -= 1
                if leading_base_count == len(ref): # ref is a prefix of alt
                    ref_border_pos = pos + len(ref) - 1
                    ref_border_seq = ref_record.seq[0:ref_border_pos]
                    alt_border_seq = ref_border_seq[0:ref_border_pos] + alt[len(ref):]
                    microhomology_size = 0
                    j = 1
                    while j < 100:
                        if ref_border_seq[-j] == alt_border_seq[-j]:
                            microhomology_size += 1
                            j += 1
                        else:
                            break
                    if microhomology_size == 0 or microhomology_size == 1:
                        df_final_variants.loc[i, "InDel_group"] = "NHEJ"
                    elif microhomology_size > 1:
                        df_final_variants.loc[i, "InDel_group"] = "MMEJ"
                    indel_size = len(alt) - len(ref) # positive for insertion, negative for deletion
                    df_final_variants.loc[i, "summary_anno"] = "MH_size=" + str(microhomology_size) + ";indel_size=" + str(indel_size)
                    df_final_variants.loc[i, "annotation"] = "Yes"
                elif leading_base_count == len(alt): # alt is a prefix of ref
                    ref_border_seq = ref_record.seq[0:(pos-1)] + ref
                    alt_border_seq = ref_border_seq[0:(pos-1)] + alt
                    microhomology_size = 0
                    j = 1
                    while j < 100:
                        if ref_border_seq[-j] == alt_border_seq[-j]:
                            microhomology_size += 1
                            j += 1
                        else:
                            break
                    if microhomology_size == 0 or microhomology_size == 1:
                        df_final_variants.loc[i, "InDel_group"] = "NHEJ"
                    elif microhomology_size > 1:
                        df_final_variants.loc[i, "InDel_group"] = "MMEJ"
                    indel_size = len(alt) - len(ref) # positive for insertion, negative for deletion
                    df_final_variants.loc[i, "summary_anno"] = "MH_size=" + str(microhomology_size) + ";indel_size=" + str(indel_size)
                    df_final_variants.loc[i, "annotation"] = "Yes"
    df_final_variants.to_excel(output_excel, index=False)


def combine_variant_anno(input_excel, output_excel):
    # combine lines with the same POS, REF, and variant_type (SNV, InDel_homopolymer, ...)
    df_final_variants = pd.read_excel(input_excel)
    output_columns = ["pos", "ref", "alt", "type", "total_count", "combined_info", "multi-allelic", "method", "sample_info", "counts", "ID_list", "fixed_REF"]
    if len(df_final_variants) == 0:
        pd.DataFrame(columns=output_columns).to_excel(output_excel, index=False)
        return
    # current columns: POS, REF, ALT, TYPE, ID_list, counts, annotation, trinucleotide_context, SNV_group, InDel_group, summary_anno
    # unused column: annotation
    df_final_variants = df_final_variants.drop(columns=["annotation"])
    combined_rows = []
    pos_list = list(df_final_variants["POS"])
    pos_list = list(set(pos_list))
    pos_list.sort()
    for pos in pos_list:
        df_tmp = df_final_variants[df_final_variants["POS"] == pos]
        df_tmp = df_tmp.reset_index(drop=True)
        SNV_rows = []
        InDel_homopolymer_rows = []
        InDel_homodimer_rows = []
        InDel_tandem_rows = []
        InDel_others_rows = []
        for i in range(len(df_tmp)):
            ref = df_tmp.loc[i, "REF"]
            alt = df_tmp.loc[i, "ALT"]
            aln_type = df_tmp.loc[i, "TYPE"]
            trinucleotide_context = df_tmp.loc[i, "trinucleotide_context"]
            SNV_group = df_tmp.loc[i, "SNV_group"]
            InDel_group = df_tmp.loc[i, "InDel_group"]
            summary_anno = df_tmp.loc[i, "summary_anno"]
            sample_info = ref + "," + alt + "," + aln_type + "," + trinucleotide_context + "," + SNV_group + "," + InDel_group + "|" + summary_anno
            ID_list = df_tmp.loc[i, "ID_list"]
            counts = df_tmp.loc[i, "counts"]
            # separate SNVs, InDels(homopolymer, homodimer, tandem), MMEJ, NHEJ, No
            if aln_type == "SNV":
                SNV_rows.append([pos, ref, "SNV", alt, sample_info, counts, ID_list])
            elif aln_type == "InDel":
                if InDel_group == "homopolymer":
                    InDel_homopolymer_rows.append([pos, ref, "InDel,homopolymer", alt, sample_info, counts, ID_list])
                elif InDel_group == "homodimer":
                    InDel_homodimer_rows.append([pos, ref, "InDel,homodimer", alt, sample_info, counts, ID_list])
                elif InDel_group == "tandem":
                    InDel_tandem_rows.append([pos, ref, "InDel,tandem", alt, sample_info, counts, ID_list])
                else:
                    InDel_others_rows.append([pos, ref, "InDel," + InDel_group, alt, sample_info, counts, ID_list])
        # combine multi-allelic SNVs: same pos, ref, different alt; lable the counts
        if len(SNV_rows) > 0:
            df_tmp_SNV = pd.DataFrame(SNV_rows)
            df_tmp_SNV.columns = ["pos", "ref", "type", "alt", "sample_info", "counts", "ID_list"]
            alt_list = list(df_tmp_SNV["alt"])
            sample_info_list = list(df_tmp_SNV["sample_info"])
            count_list = list(df_tmp_SNV["counts"])
            for i in range(len(count_list)):
                count_list[i] = str(count_list[i])
            ID_names_list = list(df_tmp_SNV["ID_list"])
            combined_rows.append([pos, ref, "SNV", "#".join(alt_list), "#".join(sample_info_list), "#".join(count_list), "#".join(ID_names_list), "code"])
        # combine multi-allelic InDels: homopolymer, homodimer, tandem: same pos, ref, repeating unit, different alt; lable the counts
        if len(InDel_homopolymer_rows) > 0:
            df_tmp_InDel_homopolymer = pd.DataFrame(InDel_homopolymer_rows)
            df_tmp_InDel_homopolymer.columns = ["pos", "ref", "type", "alt", "sample_info", "counts", "ID_list"]
            alt_list = list(df_tmp_InDel_homopolymer["alt"])
            sample_info_list = list(df_tmp_InDel_homopolymer["sample_info"])
            count_list = list(df_tmp_InDel_homopolymer["counts"])
            for i in range(len(count_list)):
                count_list[i] = str(count_list[i])
            ID_names_list = list(df_tmp_InDel_homopolymer["ID_list"])
            combined_rows.append([pos, ref, "InDel,homopolymer", "#".join(alt_list), "#".join(sample_info_list), "#".join(count_list), "#".join(ID_names_list), "code"])
        if len(InDel_homodimer_rows) > 0:
            df_tmp_InDel_homodimer = pd.DataFrame(InDel_homodimer_rows)
            df_tmp_InDel_homodimer.columns = ["pos", "ref", "type", "alt", "sample_info", "counts", "ID_list"]
            alt_list = list(df_tmp_InDel_homodimer["alt"])
            sample_info_list = list(df_tmp_InDel_homodimer["sample_info"])
            count_list = list(df_tmp_InDel_homodimer["counts"])
            for i in range(len(count_list)):
                count_list[i] = str(count_list[i])
            ID_names_list = list(df_tmp_InDel_homodimer["ID_list"])
            combined_rows.append([pos, ref, "InDel,homodimer", "#".join(alt_list), "#".join(sample_info_list), "#".join(count_list), "#".join(ID_names_list), "code"])
        if len(InDel_tandem_rows) > 0:
            df_tmp_InDel_tandem = pd.DataFrame(InDel_tandem_rows)
            df_tmp_InDel_tandem.columns = ["pos", "ref", "type", "alt", "sample_info", "counts", "ID_list"]
            alt_list = list(df_tmp_InDel_tandem["alt"])
            sample_info_list = list(df_tmp_InDel_tandem["sample_info"])
            count_list = list(df_tmp_InDel_tandem["counts"])
            for i in range(len(count_list)):
                count_list[i] = str(count_list[i])
            ID_names_list = list(df_tmp_InDel_tandem["ID_list"])
            combined_rows.append([pos, ref, "InDel,tandem", "#".join(alt_list), "#".join(sample_info_list), "#".join(count_list), "#".join(ID_names_list), "code"])
            # do not combine complex cases: MMEJ, NHEJ, No
        if len(InDel_others_rows) > 0:
            for i in range(len(InDel_others_rows)):
                pos, ref, aln_type, alt, sample_info, counts, ID_list = InDel_others_rows[i]
                combined_rows.append([pos, ref, aln_type, alt, sample_info, str(counts), ID_list, "code"])
        
    # save two copies: add a column for method: code, manual
    df_variants_combined = pd.DataFrame(combined_rows)
    df_variants_combined.columns = ["pos", "ref", "type", "alt", "sample_info", "counts", "ID_list", "method"]
    df_variants_combined = df_variants_combined.sort_values(by=["pos"])
    df_variants_combined = df_variants_combined.reset_index(drop=True)

    # add multi-allelic if there are multiple alts
    for i in range(len(df_variants_combined)):
        variant_type = df_variants_combined.loc[i, "type"]
        count = df_variants_combined.loc[i, "counts"]
        total_count = 0
        for c in count.split("#"):
            total_count += int(c)
        if "#" in count:
            df_variants_combined.loc[i, "multi-allelic"] = "multi-allelic"
        else:
            df_variants_combined.loc[i, "multi-allelic"] = "di-allelic"
        df_variants_combined.loc[i, "total_count"] = total_count
        if variant_type == "SNV":
            sample_info = df_variants_combined.loc[i, "sample_info"]
            if "#" in sample_info:
                SNV_group_list = []
                trinucleotide_context_group = ""
                for j in range(len(sample_info.split("#"))):
                    info = sample_info.split("#")[j] # # A,T,SNV,simple,AAT,A>T,-|-#A,C,SNV,simple,AAT,A>C,-|-
                    c = count.split("#")[j]
                    _, _, aln_type, trinucleotide_context, SNV_group, _ = info.split(",")
                    SNV_group_list = SNV_group_list + [SNV_group + "=" + str(c)]
                    trinucleotide_context_group = trinucleotide_context
                df_variants_combined.loc[i, "combined_info"] = trinucleotide_context_group + "," + ",".join(SNV_group_list)
            else:
                _, _, aln_type, trinucleotide_context, SNV_group, _ = sample_info.split(",")
                df_variants_combined.loc[i, "combined_info"] = trinucleotide_context + "," + SNV_group + "=" + str(total_count)
        if variant_type == "InDel,homopolymer" or variant_type == "InDel,homodimer" or variant_type == "InDel,tandem":
            sample_info = df_variants_combined.loc[i, "sample_info"] # CTTTTTTTTTT,CTTTTTTTTT,InDel,simple,-,-,homopolymer|poly-T;ref_size=10;alt_size=9
            if "#" in sample_info:
                ref_size_str_group = ""
                alt_count_list = []
                for j in range(len(sample_info.split("#"))):
                    info_1, info_2 = sample_info.split("#")[j].split("|")
                    subtype, ref_size_str, alt_size_str = info_2.split(";")
                    ref_size = ref_size_str.split("=")[1]
                    alt_size = alt_size_str.split("=")[1]
                    change_size = int(alt_size) - int(ref_size)
                    c = count.split("#")[j]
                    ref_size_str_group = ref_size_str
                    alt_count_list = alt_count_list + [str(change_size) + ":" + str(c)]
                df_variants_combined.loc[i, "combined_info"] = subtype + ";" + ref_size_str_group + ";" + ",".join(alt_count_list)
            else:
                info_1, info_2 = sample_info.split("|")
                # poly-T;ref_size=10;alt_size=9
                subtype, ref_size_str, alt_size_str = info_2.split(";")
                ref_size = ref_size_str.split("=")[1]
                alt_size = alt_size_str.split("=")[1]
                change_size = int(alt_size) - int(ref_size)
                df_variants_combined.loc[i, "combined_info"] = subtype + ";" + ref_size_str + ";" + str(change_size) + ":" + str(total_count)
        elif variant_type == "InDel,MMEJ" or variant_type == "InDel,NHEJ":
            # TTGATAATGAT,TTGAT,InDel,simple,-,-,MMEJ|MH_size=4;indel_size=-6
            sample_info = df_variants_combined.loc[i, "sample_info"]
            info_1, info_2 = sample_info.split("|")
            df_variants_combined.loc[i, "combined_info"] = info_2
        else:
            if variant_type != "SNV":
                df_variants_combined.loc[i, "combined_info"] = "-"

    # sort the columns
    df_variants_combined = df_variants_combined[["pos", "ref", "alt", "type", "total_count", "combined_info", "multi-allelic", "method", "sample_info", "counts", "ID_list"]]
    # fix the REF
    for i in range(len(df_variants_combined)):
        ref_old = df_variants_combined.loc[i, "ref"]
        multi_allelic = df_variants_combined.loc[i, "multi-allelic"]
        if multi_allelic == "di-allelic":
            sample_info = df_variants_combined.loc[i, "sample_info"]
            REF_new = sample_info.split(",")[0]
            if REF_new != ref_old:
                df_variants_combined.loc[i, "ref"] = REF_new
                df_variants_combined.loc[i, "fixed_REF"] = "Yes"
            else:
                df_variants_combined.loc[i, "fixed_REF"] = "No"
        else:
            # ACC,ACCC,InDel,-,-,homopolymer|poly-C;ref_size=2;alt_size=3#ACC,AC,InDel,-,-,homopolymer|poly-C;ref_size=2;alt_size=1
            sample_info_list = df_variants_combined.loc[i, "sample_info"].split("#")
            ref_new_list = []
            for sample_info in sample_info_list:
                ref_new_list.append(sample_info.split(",")[0])
            if len(set(ref_new_list)) == 1:
                ref_new = ref_new_list[0]
                if ref_new != ref_old:
                    df_variants_combined.loc[i, "ref"] = ref_new
                    df_variants_combined.loc[i, "fixed_REF"] = "Yes"
                else:
                    df_variants_combined.loc[i, "fixed_REF"] = "No"
            elif len(set(ref_new_list)) > 1:
                df_variants_combined.loc[i, "fixed_REF"] = "multi_ref"
                df_variants_combined.loc[i, "fixed_REF"] = "No"
    df_variants_combined.to_excel(output_excel, index=False)


def add_depth_and_frq(input_excel, variant_cov_file, output_excel, filt_excel, engine="openpxl"): # or "calamine"
    if not os.path.exists(input_excel):
        print(f"Input file {input_excel} does not exist.")
        return
    if not os.path.exists(variant_cov_file):
        print(f"Variant coverage file {variant_cov_file} does not exist.")
        return
    df = pl.read_excel(input_excel, engine=engine) # type?
    df_cov = pl.read_csv(variant_cov_file, separator="\t", has_header=False)
    df_cov = df_cov.rename({"column_1": "pos", "column_2": "depth"})
    if df.height == 0:
        df_added = df.with_columns([
            pl.lit(None, dtype=pl.Int64).alias("depth"),
            pl.lit(None, dtype=pl.Float64).alias("frequency"),
        ])
        df_added.write_excel(output_excel)
        df_added.write_excel(filt_excel)
        return
    df_added = df.join(df_cov, on="pos", how="left")
    df_added = df_added.with_columns(
        (pl.col("counts")
           .str.split("#")
           .list.eval(pl.element().cast(pl.Int64))
           .list.sum()
           / pl.col("depth"))
        .alias("frequency")
    )
    df_added.write_excel(output_excel)
    df_added_filter = df_added.filter(pl.col("frequency") >= 0.5)
    df_added_filter.write_excel(filt_excel)
    return


def run_transfer_blastn(chr_absolute_path, mito_absolute_path, plastid_absolute_path, soft_paths_dict):
    tmp_blastn_dir = "tmp_blastn_results"
    if os.path.exists(tmp_blastn_dir):
        hfbase.get_cli_output_lines("rm -rf " + tmp_blastn_dir, side_effect = True)
    os.makedirs(tmp_blastn_dir)
    # prepare mito and plastid fasta files
    hfref.replace_fasta_id("mito", mito_absolute_path, "mito.fasta")
    hfref.replace_fasta_id("plastid", plastid_absolute_path, "plastid.fasta")
    num_mito = 0
    for record in SeqIO.parse("mito.fasta", "fasta"):
        num_mito += 1
        SeqIO.write(record, "mito_" + str(num_mito) + ".fasta", "fasta")
    num_plastid = 0
    for record in SeqIO.parse("plastid.fasta", "fasta"):
        num_plastid += 1
        SeqIO.write(record, "plastid_" + str(num_plastid) + ".fasta", "fasta")
    # prepare chr fasta file
    hfref.replace_fasta_id("Chr", chr_absolute_path, "Chr.fasta")
    Chr_records = list(SeqIO.parse("Chr.fasta", "fasta"))
    num_Chr = len(Chr_records)
    for i in range(num_mito):
        genome_absolute_path = "mito_" + str(i+1) + ".fasta"
        results_filename = "mito_" + str(i+1) + "_blastn_results.txt"
        if os.path.exists(tmp_blastn_dir + "/" + results_filename):
            os.remove(tmp_blastn_dir + "/" + results_filename)
        with cf.ThreadPoolExecutor(num_Chr) as tex:
            futures = [tex.submit(hfvar.run_blastn_sorter_single, Chr_record, genome_absolute_path, results_filename, soft_paths_dict, tmp_blastn_dir) for Chr_record in Chr_records]
            results = [future.result() for future in cf.as_completed(futures)]
        # convert to excel
        blastn_lines = hfbase.get_file_lines(tmp_blastn_dir + "/" + results_filename)
        for blastn_line in blastn_lines:
            df_blastn_alignments = pd.DataFrame(columns=["chr_name", "chr_len", "aln_index", "aln_len", "aln_olp_len", "aln_idt", "aln_strand", "aln_qs", "aln_qe", "aln_ss", "aln_se", "aln_cn"])
            chr_name, _, chr_len, aln_type, percent_total, *align_info_list_of_dict = blastn_line.split("\t")
            for j in range(len(align_info_list_of_dict)):
                align_info_list_of_dict[j] = OrderedDict([ (x.split("=")[0], x.split("=")[1]) for x in align_info_list_of_dict[j].split(";") ])
                aln_index = align_info_list_of_dict[j]["aln"] # a string
                aln_len = align_info_list_of_dict[j]["len"]
                aln_olp_len = align_info_list_of_dict[j]["olp"]
                aln_idt = align_info_list_of_dict[j]["idt"]
                aln_strand = align_info_list_of_dict[j]["strand"]
                aln_qs = align_info_list_of_dict[j]["qs"]
                aln_qe = align_info_list_of_dict[j]["qe"]
                aln_ss = align_info_list_of_dict[j]["ss"]
                aln_se = align_info_list_of_dict[j]["se"]
                aln_cn = align_info_list_of_dict[j]["cn"]
                df_blastn_alignments.loc[j] = [chr_name, chr_len, aln_index, aln_len, aln_olp_len, aln_idt, aln_strand, aln_qs, aln_qe, aln_ss, aln_se, aln_cn]
            df_blastn_alignments.to_excel(chr_name + "_vs_" + "mito_" + str(i+1) + "_blastn_alignments.xlsx", index=False)

    for i in range(num_plastid):
        genome_absolute_path = "plastid_" + str(i+1) + ".fasta"
        results_filename = "plastid_" + str(i+1) + "_blastn_results.txt"
        if os.path.exists(tmp_blastn_dir + "/" + results_filename):
            os.remove(tmp_blastn_dir + "/" + results_filename)
        with cf.ThreadPoolExecutor(num_Chr) as tex:
            futures = [tex.submit(hfvar.run_blastn_sorter_single, Chr_record, genome_absolute_path, results_filename, soft_paths_dict, tmp_blastn_dir) for Chr_record in Chr_records]
            results = [future.result() for future in cf.as_completed(futures)]
        # convert to excel
        blastn_lines = hfbase.get_file_lines(tmp_blastn_dir + "/" + results_filename)
        for blastn_line in blastn_lines:
            df_blastn_alignments = pd.DataFrame(columns=["chr_name", "chr_len", "aln_index", "aln_len", "aln_olp_len", "aln_idt", "aln_strand", "aln_qs", "aln_qe", "aln_ss", "aln_se", "aln_cn"])
            chr_name, _, chr_len, aln_type, percent_total, *align_info_list_of_dict = blastn_line.split("\t")
            for j in range(len(align_info_list_of_dict)):
                align_info_list_of_dict[j] = OrderedDict([ (x.split("=")[0], x.split("=")[1]) for x in align_info_list_of_dict[j].split(";") ])
                aln_index = align_info_list_of_dict[j]["aln"] # a string
                aln_len = align_info_list_of_dict[j]["len"]
                aln_olp_len = align_info_list_of_dict[j]["olp"]
                aln_idt = align_info_list_of_dict[j]["idt"]
                aln_strand = align_info_list_of_dict[j]["strand"]
                aln_qs = align_info_list_of_dict[j]["qs"]
                aln_qe = align_info_list_of_dict[j]["qe"]
                aln_ss = align_info_list_of_dict[j]["ss"]
                aln_se = align_info_list_of_dict[j]["se"]
                aln_cn = align_info_list_of_dict[j]["cn"]
                df_blastn_alignments.loc[j] = [chr_name, chr_len, aln_index, aln_len, aln_olp_len, aln_idt, aln_strand, aln_qs, aln_qe, aln_ss, aln_se, aln_cn]
            df_blastn_alignments.to_excel(chr_name + "_vs_" + "plastid_" + str(i+1) + "_blastn_alignments.xlsx", index=False)
    hfbase.get_cli_output_lines("rm -rf " + tmp_blastn_dir, side_effect = True)
    return num_Chr, num_mito, num_plastid


def merge_fragments_files(chr_name, num_mito, num_plastid):
    mito_fragments_df_list = []
    for j in range(num_mito):
        mito_contig_file = chr_name + "_vs_mito_" + str(j + 1) + "_blastn_alignments.xlsx"
        if os.path.exists(mito_contig_file):
            df_mito = pl.read_excel(mito_contig_file, engine="calamine")
            df_mito = df_mito.with_columns(
                pl.lit("mito_" + str(j + 1)).alias("contig_name")
            )
            df_mito = df_mito.select(
                pl.col("chr_name"),
                pl.col("contig_name"),
                pl.col("aln_len"),
                pl.col("aln_idt"),
                pl.col("aln_strand"),
                pl.col("aln_qs"),
                pl.col("aln_qe"),
                pl.col("aln_ss"),
                pl.col("aln_se")
            )
            mito_fragments_df_list.append(df_mito)
    try:
        mito_fragments_df = pl.concat(mito_fragments_df_list)
    except ValueError:
        mito_fragments_df = pl.DataFrame(
            {
                "chr_name": [],
                "contig_name": [],
                "aln_len": [],
                "aln_idt": [],
                "aln_strand": [],
                "aln_qs": [],
                "aln_qe": [],
                "aln_ss": [],
                "aln_se": []
            }
        )
    plastid_fragments_df_list = []
    for k in range(num_plastid):
        plastid_contig_file = chr_name + "_vs_plastid_" + str(k + 1) + "_blastn_alignments.xlsx"
        if os.path.exists(plastid_contig_file):
            df_plastid = pl.read_excel(plastid_contig_file, engine="calamine")
            df_plastid = df_plastid.with_columns(
                pl.lit("plastid_" + str(k + 1)).alias("contig_name")
            )
            df_plastid = df_plastid.select(
                pl.col("chr_name"),
                pl.col("contig_name"),
                pl.col("aln_len"),
                pl.col("aln_idt"),
                pl.col("aln_strand"),
                pl.col("aln_qs"),
                pl.col("aln_qe"),
                pl.col("aln_ss"),
                pl.col("aln_se")
            )
            plastid_fragments_df_list.append(df_plastid)
    try:
        plastid_fragments_df = pl.concat(plastid_fragments_df_list)
    except ValueError:
        plastid_fragments_df = pl.DataFrame(
            {
                "chr_name": [],
                "contig_name": [],
                "aln_len": [],
                "aln_idt": [],
                "aln_strand": [],
                "aln_qs": [],
                "aln_qe": [],
                "aln_ss": [],
                "aln_se": []
            }
        )
    try:
        all_fragments_df = pl.concat([mito_fragments_df, plastid_fragments_df])
    except ValueError:
        all_fragments_df = pl.DataFrame(
            {
                "chr_name": [],
                "contig_name": [],
                "aln_len": [],
                "aln_idt": [],
                "aln_strand": [],
                "aln_qs": [],
                "aln_qe": [],
                "aln_ss": [],
                "aln_se": []
            }
        )
    all_fragments_df = all_fragments_df.with_columns(
        pl.col("chr_name").cast(pl.Utf8),
        pl.col("contig_name").cast(pl.Utf8),
        pl.col("aln_len").cast(pl.Int64),
        pl.col("aln_idt").cast(pl.Float64),
        pl.col("aln_strand").cast(pl.Utf8),
        pl.col("aln_qs").cast(pl.Int64),
        pl.col("aln_qe").cast(pl.Int64),
        pl.col("aln_ss").cast(pl.Int64),
        pl.col("aln_se").cast(pl.Int64)
    )
    return all_fragments_df


def remove_inner_fragments(all_fragments_df):
    num_removed = 0
    for i in range(all_fragments_df.height):
        if all_fragments_df[i, "keep_or_remove"] == "remove":
            continue
        for j in range(i + 1, all_fragments_df.height):
            if all_fragments_df[j, "keep_or_remove"] == "remove":
                continue
            if all_fragments_df[i, "aln_qs"] == all_fragments_df[j, "aln_qs"] and all_fragments_df[i, "aln_qe"] == all_fragments_df[j, "aln_qe"]:
                num_removed += 1
                all_fragments_df[i, "has_inner_fragment"] = True
                if all_fragments_df[i, "inner_fragments"] == "":
                    all_fragments_df[i, "inner_fragments"] = all_fragments_df[j, "contig_name"] + "=" + str(all_fragments_df[j, "aln_qs"]) + "," + str(all_fragments_df[j, "aln_qe"])
                else:
                    all_fragments_df[i, "inner_fragments"] = ",".join(all_fragments_df[i, "inner_fragments"].split(",") + [all_fragments_df[j, "contig_name"] + "=" + str(all_fragments_df[j, "aln_qs"]) + "," + str(all_fragments_df[j, "aln_qe"])])
                all_fragments_df[i, "keep_or_remove"] = "same"
                all_fragments_df[j, "keep_or_remove"] = "remove"
                return num_removed, all_fragments_df
            else:
                if all_fragments_df[i, "aln_qs"] <= all_fragments_df[j, "aln_qs"] and all_fragments_df[i, "aln_qe"] >= all_fragments_df[j, "aln_qe"]:
                    num_removed += 1
                    all_fragments_df[i, "has_inner_fragment"] = True
                    if all_fragments_df[i, "inner_fragments"] == "":
                        all_fragments_df[i, "inner_fragments"] = all_fragments_df[j, "contig_name"] + "=" + str(all_fragments_df[j, "aln_qs"]) + "," + str(all_fragments_df[j, "aln_qe"])
                    else:
                        all_fragments_df[i, "inner_fragments"] = ",".join(all_fragments_df[i, "inner_fragments"].split(",") + [all_fragments_df[j, "contig_name"] + "=" + str(all_fragments_df[j, "aln_qs"]) + "," + str(all_fragments_df[j, "aln_qe"])])
                    all_fragments_df[i, "keep_or_remove"] = "keep"
                    all_fragments_df[j, "keep_or_remove"] = "remove"
                    return num_removed, all_fragments_df
    return num_removed, all_fragments_df


def merge_numt_nupt(num_Chr, num_mito, num_plastid):
    for i in range(num_Chr):
        chr_name = "Chr_" + str(i + 1)
        all_fragments_df = merge_fragments_files(chr_name, num_mito, num_plastid) # chr_name, contig_name, aln_len, aln_idt, aln_strand, aln_qs, aln_qe, aln_ss, aln_se
        all_fragments_df = all_fragments_df.sort("aln_qs", "contig_name", descending=[False, False])
        all_fragments_df = all_fragments_df.with_columns(
            (pl.arange(1, all_fragments_df.height + 1)).alias("index").cast(pl.Int64), # to record the original index of the row
            pl.lit(False).alias("has_inner_fragment").cast(pl.Boolean), # to record whether the fragment has inner fragments
            pl.lit("keep").alias("keep_or_remove").cast(pl.Utf8), # to record whether the fragment should be kept or removed
            pl.lit("").alias("inner_fragments").cast(pl.Utf8) # to record the inner fragments
        )
        num_removed, all_fragments_df = remove_inner_fragments(all_fragments_df)
        while num_removed > 0:
            num_removed, all_fragments_df = remove_inner_fragments(all_fragments_df)
        all_fragments_df = all_fragments_df.filter((pl.col("keep_or_remove") == "keep") | (pl.col("keep_or_remove") == "same"))
        all_fragments_df = all_fragments_df.with_columns(
            pl.lit(0).alias("aln_olp_len").cast(pl.Int64), # to record the overlap length
        )
        for i in range(all_fragments_df.height):
            if i != (all_fragments_df.height - 1):
                all_fragments_df[i, "aln_olp_len"] = all_fragments_df[i, "aln_qe"] - all_fragments_df[i + 1, "aln_qs"] + 1
        all_fragments_df = all_fragments_df.with_columns(
            pl.col("index").alias("original_index"),
            (pl.arange(1, all_fragments_df.height + 1)).alias("index").cast(pl.Int64),
        )
        all_fragments_df.write_excel(chr_name + "_all_fragments.xlsx")


def get_raw_transfer_groups(num_Chr):
    for i in range(num_Chr):
        chr_name = "Chr_" + str(i + 1)
        df_blastn_alignments = pl.read_excel(chr_name + "_all_fragments.xlsx", engine="calamine")
        # if df_blastn_alignments is empty, skip
        if df_blastn_alignments.height == 0:
            continue
        index = 1.000
        threshold = -10000 # set or change the threshold for grouping
        df_blastn_alignments = df_blastn_alignments.with_columns(
            pl.col("index").cast(pl.Float64)
        )
        for j in range(df_blastn_alignments.height):
            if j != df_blastn_alignments.height - 1:
                if df_blastn_alignments[j, "aln_olp_len"] < threshold:
                    df_blastn_alignments[j, "index"] = index
                    index = int(index) + 1
                elif df_blastn_alignments[j, "aln_olp_len"] >= threshold:
                    df_blastn_alignments[j, "index"] = index
                    index = index + 0.001
            else:
                df_blastn_alignments["index"][j] = index # the last fragment
        df_blastn_alignments = df_blastn_alignments.with_columns(
            group_index = pl.col("index").cast(pl.Int64)
        )
        # use pandas to calculate group features
        df_blastn_alignments = df_blastn_alignments.to_pandas()
        # calculate group length for each transfer event
        groups = list(set([ int(group) for group in df_blastn_alignments["group_index"] ]))
        group_info_dict = dict()
        for group in groups:
            df_group = df_blastn_alignments[df_blastn_alignments["group_index"] == group]
            group_start = min(df_group["aln_qs"])
            group_end = max(df_group["aln_qe"])
            group_size = group_end - group_start + 1
            group_info_dict[group] = [group_start, group_end, group_size]
        for j in df_blastn_alignments.index:
            df_blastn_alignments.loc[j, "group_start"] = group_info_dict[df_blastn_alignments.loc[j, "group_index"]][0]
            df_blastn_alignments.loc[j, "group_end"] = group_info_dict[df_blastn_alignments.loc[j, "group_index"]][1]
            df_blastn_alignments.loc[j, "group_size"] = group_info_dict[df_blastn_alignments.loc[j, "group_index"]][2]
            if df_blastn_alignments.loc[j, "group_size"] < 100:
                df_blastn_alignments.loc[j, "group_size_flag"] = "small"
            elif df_blastn_alignments.loc[j, "group_size"] >= 100 and df_blastn_alignments.loc[j, "group_size"] < 1000:
                df_blastn_alignments.loc[j, "group_size_flag"] = "medium"
            elif df_blastn_alignments.loc[j, "group_size"] >= 1000 and df_blastn_alignments.loc[j, "group_size"] < 10000:
                df_blastn_alignments.loc[j, "group_size_flag"] = "large"
            elif df_blastn_alignments.loc[j, "group_size"] >= 10000 and df_blastn_alignments.loc[j, "group_size"] < 100000:
                df_blastn_alignments.loc[j, "group_size_flag"] = "huge"
            elif df_blastn_alignments.loc[j, "group_size"] >= 100000:
                df_blastn_alignments.loc[j, "group_size_flag"] = "giant"
        # determine_max_fragment_size
        group_index_list = list(set(list(df_blastn_alignments["group_index"])))
        group_index_list.sort()
        for group_index in group_index_list:
            df_group = df_blastn_alignments[df_blastn_alignments["group_index"] == group_index]
            max_fragment_size = max(list(df_group["aln_len"]))
            for i in range(len(df_group)):
                df_blastn_alignments.loc[df_group.index[i], "max_fragment_size"] = max_fragment_size
        # fragments in the same group overlap or less than 100-bp away are chimeric
        for j in range(len(df_blastn_alignments)):
            if j == (len(df_blastn_alignments) - 1) and j != 0:
                group_index_j = df_blastn_alignments.loc[j, "group_index"]
                group_index_prev = df_blastn_alignments.loc[j-1, "group_index"]
                if group_index_j == group_index_prev:
                    df_blastn_alignments.loc[j, "chimeric"] = df_blastn_alignments.loc[j-1, "chimeric"]
            else:
                group_index_j = df_blastn_alignments.loc[j, "group_index"]
                group_index_next = df_blastn_alignments.loc[j+1, "group_index"]
                aln_olp_len_j = df_blastn_alignments.loc[j, "aln_olp_len"]
                if group_index_j == group_index_next and abs(aln_olp_len_j) <= 100:
                    df_blastn_alignments.loc[j, "chimeric"] = "yes"
                else:
                    df_blastn_alignments.loc[j, "chimeric"] = "no"
        # calculate nuclear percentage: sum of all aln_olp_len with positive values / group_size
        for j in range(len(df_blastn_alignments)):
            group_index = df_blastn_alignments.loc[j, "group_index"]
            df_group = df_blastn_alignments[df_blastn_alignments["group_index"] == group_index]
            sum_aln_olp_len = 0
            for k in df_group.index:
                if df_group.loc[k, "aln_olp_len"] < 0:
                    sum_aln_olp_len += abs(df_group.loc[k, "aln_olp_len"])
            # remove the last negative overlap length
            if df_group.iloc[-1]["aln_olp_len"] < 0:
                sum_aln_olp_len -= abs(df_group.iloc[-1]["aln_olp_len"])
            group_size = df_blastn_alignments.loc[j, "group_size"]
            nuclear_percentage = sum_aln_olp_len / group_size
            df_blastn_alignments.loc[j, "N_percentage"] = nuclear_percentage
        # determine group type depending on "contig_name", "keep_or_remove"
        # first mito, plastid, or mixed; then same
        for j in range(len(df_blastn_alignments)):
            group_index = df_blastn_alignments.loc[j, "group_index"]
            df_group = df_blastn_alignments[df_blastn_alignments["group_index"] == group_index]
            contig_name_set = set([item.split("_")[0] for item in list(df_group["contig_name"])])
            if "mito" in contig_name_set and "plastid" not in contig_name_set:
                df_blastn_alignments.loc[j, "group_type"] = "mito"
            elif "plastid" in contig_name_set and "mito" not in contig_name_set:
                df_blastn_alignments.loc[j, "group_type"] = "plastid"
            elif "mito" in contig_name_set and "plastid" in contig_name_set:
                df_blastn_alignments.loc[j, "group_type"] = "mixed"
            if "same" in list(df_group["keep_or_remove"]):
                df_blastn_alignments.loc[j, "group_type"] = "same"
        df_blastn_alignments = df_blastn_alignments[["index", "group_index", "group_type", "group_start", "group_end", "group_size", "group_size_flag", "max_fragment_size", "chimeric", "N_percentage", "chr_name", "contig_name", "aln_len", "aln_olp_len", "aln_idt", "aln_strand", "aln_qs", "aln_qe", "aln_ss", "aln_se", "has_inner_fragment", "keep_or_remove", "inner_fragments"]]
        df_blastn_alignments.to_excel("grouped_" + chr_name + "_blastn_alignments.xlsx", index=False)
