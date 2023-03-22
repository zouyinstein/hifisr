import pandas as pd
from Bio import SeqIO
import sys
import os

if not os.path.exists(sys.argv[1]):
    print("No variants. Use the old reference.")
    exit()

df_bcftools = pd.read_table(sys.argv[1], header=None)  # AH-7_mito_bcftools_checked.txt
old_ref = SeqIO.read(sys.argv[2], "fasta") # AH-7_mito_v0.fasta
new_seq = ""
start = 0
POS_REF_ALT_list = [ (df_bcftools.loc[i, 0], df_bcftools.loc[i, 1], df_bcftools.loc[i, 2].split(",")[0]) for i in range(len(df_bcftools)) if df_bcftools.loc[i, 3].split(":")[0] == "1/1" ]
# [(120654, 'A', 'AGTGGC'), (243635, 'ACC', 'ACCC'), (250222, 'AG', 'AGG'), (308606, 'AA', 'AAGCCA')]
for i in range(len(POS_REF_ALT_list)):
    new_seq += old_ref.seq[start:(POS_REF_ALT_list[i][0]-1)]
    new_seq += POS_REF_ALT_list[i][2]
    print("replace " + POS_REF_ALT_list[i][1] + " with " + POS_REF_ALT_list[i][2] + " at position " + str(POS_REF_ALT_list[i][0]))
    start = POS_REF_ALT_list[i][0] + len(POS_REF_ALT_list[i][1]) - 1
new_seq += old_ref.seq[start:len(old_ref.seq)]

with open("replaced_by_bcftools_" + sys.argv[2], "wt") as fout:
    print(">" + old_ref.id + "_replaced_by_bcftools", file=fout)
    print(new_seq, file=fout)
print("old ref length: " + str(len(old_ref.seq)))
print("new ref length: " + str(len(new_seq)))