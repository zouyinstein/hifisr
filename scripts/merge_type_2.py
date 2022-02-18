import pandas as pd
import re
import sys
import os
import subprocess


sub_type = sys.argv[1]
samples_txt_prefix = sys.argv[2]
genome = sys.argv[3]
with open("../../../" + samples_txt_prefix + ".txt", "rt") as fin:
    samples = [sample.rstrip("\n") for sample in fin.readlines()]

data_count = pd.read_table("count_se1_ss2_" + sub_type + "_" + samples_txt_prefix + "_" + genome + ".txt", header=None)
all_se1_ss2 = {(data_count[0][i], data_count[1][i]) for i in range(len(data_count))}
all_list = list(all_se1_ss2)
all_list.sort(key=lambda x:x[0])
df_count=pd.DataFrame(columns=["se1", "ss2"] + samples, index=all_list)
for i in range(len(data_count)):
    index_tuple = (data_count[0][i], data_count[1][i])
    df_count[data_count[4][i]][index_tuple] = data_count[2][i]
df_count.fillna(0, inplace=True)
df_count.to_excel(samples_txt_prefix + "_se1_ss2_annotation.xlsx")
df_count = pd.read_excel(samples_txt_prefix + "_se1_ss2_annotation.xlsx")

df_count.rename(columns={"Unnamed: 0":"se1_ss2"}, inplace=True)
data1 = pd.read_table("count_se1_ss2_" + sub_type + "_" + samples_txt_prefix + "_" + genome + ".txt", header=None)
for i in range(len(df_count)):
    cords = re.findall("\d+", df_count.loc[i, "se1_ss2"])
    data1_tmp = data1[(data1[0] == int(cords[0])) & (data1[1] == int(cords[1]))]
    df_count.loc[i, "se1"] = int(cords[0])
    df_count.loc[i, "ss2"] = int(cords[1])
    df_count.loc[i, "index"] = i + 1
    df_count.loc[i, "repsize"] = str(data1_tmp.iloc[0].iat[3]).lstrip()
    df_count.loc[i, "min_repsize"] = int(re.findall("\d+", df_count.loc[i, "repsize"])[0])
    df_count.loc[i, "max_repsize"] = int(re.findall("\d+", df_count.loc[i, "repsize"])[-1])
    df_count.loc[i, "mid_repsize"] = (df_count.loc[i, "min_repsize"] + df_count.loc[i, "max_repsize"])/2
df_count = df_count[["index", "se1_ss2", "se1", "ss2", "repsize", "min_repsize", "mid_repsize", "max_repsize"] + samples]
df_count.sort_values("mid_repsize", inplace = True, ascending=False)
df_count.to_excel(samples_txt_prefix + "_se1_ss2_annotation.xlsx", index=False)
df_count.to_excel(samples_txt_prefix + "_se1_ss2_backup.xlsx", index=False)

for sample in samples:
    for i in range(len(df_count)):
        data = pd.read_table("../../../" + sample + "/" + genome + "/" + sub_type + "/" + sub_type + "_" + sample + "_" + genome + ".txt", header=None)  
        data_index = data[(data[8] == df_count.loc[i, "se1"]) & (data[12] == df_count.loc[i, "ss2"])]
        ids = [id.split(":")[1] for id in list(data_index[0])]
        if not os.path.exists("reads_group_ids/" + str(int(df_count.loc[i, "index"]))):
            ret = subprocess.call("mkdir -p " + "reads_group_ids/" + str(int(df_count.loc[i, "index"])), shell=True)
        with open("reads_group_ids/" + str(int(df_count.loc[i, "index"])) + "/" + sample + "_ids.list", "wt") as fout:
            print("\n".join(ids), file=fout)
