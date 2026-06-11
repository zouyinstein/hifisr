import pandas as pd
from Bio import SeqIO
import sys
import _bootstrap  # noqa: F401
import os


def correct_pos_dict(old_fasta, new_fasta, info_dict):
    # can deal with 2 overlapping variants, not 3 overlapping variants, or more
    record = SeqIO.read(old_fasta, "fasta")
    number_of_variants = len(info_dict)
    POS_list = list(info_dict.keys())
    POS_list.sort()
    new_fasta_seq = ""
    if number_of_variants == 1:
        POS_1 = POS_list[0]
        REF_1, ALT_1 = info_dict[POS_1]
        new_fasta_seq = record.seq[:POS_1-1] + ALT_1 + record.seq[POS_1-1+len(REF_1):]
    elif number_of_variants == 2:
        POS_1 = POS_list[0]
        REF_1, ALT_1 = info_dict[POS_1]
        POS_2 = POS_list[1]
        REF_2, ALT_2 = info_dict[POS_2]
        delta = POS_1-1+len(REF_1) - (POS_2-1)
        if delta <= 0:
            replace_seq = ALT_1 + record.seq[POS_1-1+len(REF_1):POS_2-1] + ALT_2
        elif delta > 0:
            replace_seq = ALT_1[0:(len(ALT_1)-delta)] + ALT_2
        new_fasta_seq = record.seq[:POS_1-1] + replace_seq + record.seq[POS_2-1+len(REF_2):]
    else:
        REF_list = [ info_dict[POS][0] for POS in POS_list ]
        ALT_list = [ info_dict[POS][1] for POS in POS_list ]
        new_fasta_seq = ""
        replace_seq = ""
        # POS_1-1+len(REF_1) - (POS_2-1)
        delta_list = [ POS_list[i]-1+len(REF_list[i]) - (POS_list[i+1]-1) for i in range(number_of_variants-1) ]
        for i in range(number_of_variants):
            if i == 0:
                if delta_list[i] > 0:
                    REF_list[i] = REF_list[i][:len(REF_list[i])-delta_list[i]] # remove the overlapping part
                    ALT_list[i] = ALT_list[i][:len(ALT_list[i])-delta_list[i]] # remove the overlapping part
                replace_seq = ALT_list[i] + record.seq[POS_list[i]-1+len(REF_list[i]):POS_list[i+1]-1]
                new_fasta_seq = record.seq[:POS_list[i]-1] + replace_seq
            elif i == number_of_variants - 1:
                new_fasta_seq = new_fasta_seq + ALT_list[i] + record.seq[POS_list[i]-1+len(REF_list[i]):]
            else:
                if delta_list[i] > 0:
                    REF_list[i] = REF_list[i][:len(REF_list[i])-delta_list[i]] # remove the overlapping part
                    ALT_list[i] = ALT_list[i][:len(ALT_list[i])-delta_list[i]] # remove the overlapping part
                replace_seq = ALT_list[i] + record.seq[POS_list[i]-1+len(REF_list[i]):POS_list[i+1]-1]
                new_fasta_seq = new_fasta_seq + replace_seq

    with open(new_fasta, "wt") as fout:
        print(">" + record.id + "_converted", file=fout)
        print(new_fasta_seq, file=fout)
    print("old length:", len(record.seq))
    print("new length:", len(new_fasta_seq))
    return


def get_file_lines(file):
    with open(file, "rt") as fin:
        lines = [ line.rstrip("\n") for line in fin.readlines() ]
    return lines


old_fasta = sys.argv[1]
new_fasta = sys.argv[2]
variant_file = sys.argv[3] # pos, ref, alt

pos_ref_alt_lines = get_file_lines(variant_file)
if len(pos_ref_alt_lines) == 0:
    os.system(f"cp {old_fasta} {new_fasta}")
    sys.exit(0)
info_dict = {}
for line in pos_ref_alt_lines:
    pos, ref, alt = line.split("\t")
    info_dict[int(pos)] = (ref, alt)
correct_pos_dict(old_fasta, new_fasta, info_dict)
