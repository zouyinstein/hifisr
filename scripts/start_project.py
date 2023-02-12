import pandas as pd
import sys
import os

# input_files.txt
df_samples = pd.read_table(sys.argv[1], header=None)    
threads = sys.argv[2]

with open("work.sh", "wt") as fout:
    print("ls start*.sh | while read i; do bash $i; done", file=fout)

for i in range(len(df_samples)):
    sample = df_samples.loc[i, 0]
    reads = df_samples.loc[i, 1]
    ref_fa = df_samples.loc[i, 2]
    mito_fa = df_samples.loc[i, 3]
    plastid_fa = df_samples.loc[i, 4]
    os.mkdir(sample)
    
    with open("start_" + sample + ".sh", "at") as fout:
        print("cd " + sample, file=fout)
        print("ln -sf " + reads + " " + sample + ".fastq", file=fout)  # fastq
        print("ln -sf " + ref_fa + " " + sample + "_ref.fa", file=fout) 
        print("ln -sf " + mito_fa + " " + sample + "_mito.fa", file=fout)  
        print("ln -sf " + plastid_fa + " " + sample + "_plastid.fa", file=fout)
        print("bash work.sh", file=fout)
        print("cd ..", file=fout)
   
    with open(sample + "/work.sh", "at") as fout:
        print("python ../../hifisr.py -s " + sample + " -t " + threads + " -i fastq single > $(date +%s).log 2> $(date +%s).err", file=fout)
        print("bash ../../scripts/run_mecat_flye.sh " + sample + " mito > assembly_mito_$(date +%s).log 2> assembly_mito_$(date +%s).err" + threads, file=fout)
        print("bash ../../scripts/run_mecat_flye.sh " + sample + " plastid > assembly_plastid_$(date +%s).log 2> assembly_plastid_$(date +%s).err" + threads, file=fout)

