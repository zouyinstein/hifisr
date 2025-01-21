import pandas as pd
import pysam
from Bio import SeqIO
import sys
import os

if os.path.getsize(sys.argv[1]) == 0:
    print("No variants detected by bcftools.")
    exit()

df_bcftools = pd.read_table(sys.argv[1], header=None)  # AH-7_mito_bcftools_small.vcf
samfile = pysam.AlignmentFile(sys.argv[2], "rb") # type_1_AH-7.mito_filt_v0.sorted.bam
ref_ID = SeqIO.read(sys.argv[3], "fasta").id # AH-7_mito_v0.fasta
if not samfile.check_index():
    print("No bam index")
    exit()
    
for i in range(len(df_bcftools)):
    POS = df_bcftools.loc[i, 0]
    REF = df_bcftools.loc[i, 1]
    ALT = df_bcftools.loc[i, 2]  # could be a list "A,T"
    ALT_length_list = [ len(item) for item in ALT.split(",") ]
    BCFTOOLS = df_bcftools.loc[i, 3] 
    if len(REF) == 1 and max(ALT_length_list) == 1:   # SNV
        counts_array = samfile.count_coverage(contig=ref_ID, start=POS-1, stop=POS, quality_threshold=0)    # ACGT, N
        counts_str = str(list(counts_array[0])[0]) + "," + str(list(counts_array[1])[0]) + ","+ str(list(counts_array[2])[0]) + "," + str(list(counts_array[3])[0])
        print(str(POS) + "\t" + REF + "\t"+ ALT + "\t" + BCFTOOLS + "\t" + counts_str)
    else:  # indel
        counts_dict = dict()
        for pc in samfile.pileup(contig=ref_ID, start=POS-1, stop=POS):
            if pc.pos == (POS-1):
                for pr in pc.pileups:
                    if counts_dict.get(pr.indel) == None:
                        counts_dict[pr.indel] = 1
                    else:
                        counts_dict[pr.indel] += 1
        counts_keys = list(counts_dict.keys())
        counts_keys.sort()
        counts_list_info = [ str(item) + ":" + str(counts_dict[item]) for item in counts_keys ]
        print(str(POS) + "\t" + REF + "\t"+ ALT + "\t" + BCFTOOLS + "\t" +  ";".join(counts_list_info))
samfile.close()

