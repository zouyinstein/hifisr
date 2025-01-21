import pandas as pd
import sys
import os

# input_files.txt
df_samples = pd.read_table(sys.argv[1], header=None)    
threads = sys.argv[2]

with open("work.sh", "wt") as fout:
    print("ls start*.sh | while read i; do bash $i; done", file=fout)
with open("clean.sh", "wt") as fout:
    print("cut -f1 " + sys.argv[1] + " | while read i; do rm -rf $i start_${i}.sh work.sh clean.sh; done", file=fout)

for i in range(len(df_samples)):
    sample = df_samples.loc[i, 0]
    mito_ref = df_samples.loc[i, 1]
    mito_reads = df_samples.loc[i, 2]
    plastid_ref = df_samples.loc[i, 3]
    plastid_reads = df_samples.loc[i, 4]
    if not os.path.exists(sample):
        os.mkdir(sample)
    
    with open("start_" + sample + ".sh", "at") as fout:
        print("cd " + sample, file=fout)
        print("ln -sf " + mito_ref + " mito_ref.fasta", file=fout)
        print("ln -sf " + mito_reads + " mito_reads.fastq", file=fout)
        print("ln -sf " + plastid_ref + " plastid_ref.fasta", file=fout)
        print("ln -sf " + plastid_reads + " plastid_reads.fastq", file=fout)
        print("bash work.sh", file=fout)
        print("cd ..", file=fout)
   
    with open(sample + "/work.sh", "at") as fout:
        print("bash ../run_polish_by_bcftools.sh mito_ref.fasta mito_reads.fastq " + threads + " > mito_$(date +%s).log 2> mito_$(date +%s).err", file=fout)  
        print("bash ../run_polish_by_bcftools.sh plastid_ref.fasta plastid_reads.fastq " + threads + " > plastid_$(date +%s).log 2> plastid_$(date +%s).err", file=fout) 

