import pandas as pd
import re
import sys
import os

blast_file = sys.argv[1]
sample=sys.argv[2]
data = pd.read_table(blast_file, header=None)

all_se1_ss2 = { (data[8][i], data[12][i]) for i in range(len(data)) }
if not os.path.exists("reads_group"):
    os.mkdir("reads_group")
for item in list(all_se1_ss2):
    data_tmp = data[(data[8] == item[0]) & (data[12] == item[1])]
    data_tmp.to_csv("reads_group/" + str(item[0]) + "_" + str(item[1]) + "_" + blast_file, header=False, index=False, sep="\t")
    repsize = list(set(data_tmp[3]))
    repsize.sort()
    repsize_str_list = [str(item) for item in repsize]
    repsize_str = ",".join(repsize_str_list)
    with open("count_se1_ss2_" + blast_file, "at") as f_count:
        print("\t".join([str(item[0]), str(item[1]), str(len(data_tmp)), repsize_str, sample]), file=f_count)

data_count = pd.read_table("count_se1_ss2_" + blast_file, header=None)
all_se1_ss2 = {(data_count[0][i], data_count[1][i]) for i in range(len(data_count))}
all_list = list(all_se1_ss2)
all_list.sort(key=lambda x:x[0])

df_count=pd.DataFrame(columns=["se1", "ss2", sample], index=all_list)
for i in range(len(data_count)):
    index_tuple = (data_count[0][i], data_count[1][i])
    df_count[data_count[4][i]][index_tuple] = data_count[2][i]
df_count.fillna(0, inplace=True)
df_count.to_excel(sample + "_se1_ss2_annotation.xlsx")
df_count = pd.read_excel(sample + "_se1_ss2_annotation.xlsx")
df_count.rename(columns={"Unnamed: 0":"se1_ss2"}, inplace=True)
data1 = pd.read_table("count_se1_ss2_" + blast_file, header=None)
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
df_count = df_count[["index", "se1_ss2", "se1", "ss2", "repsize", "min_repsize", "mid_repsize", "max_repsize", sample]]
df_count.sort_values("mid_repsize", inplace = True, ascending=False)
df_count.to_excel(sample + "_se1_ss2_annotation.xlsx", index=False)
df_count.to_excel(sample + "_se1_ss2_backup.xlsx", index=False)
if not os.path.exists("reads_group_ids"):
    os.mkdir("reads_group_ids")
for i in range(len(df_count)):
    data_index = data[(data[8] == df_count.loc[i, "se1"]) & (data[12] == df_count.loc[i, "ss2"])]
    with open("reads_group_ids/" + str(int(df_count.loc[i, "index"])) + "_ids.list", "wt") as fout:
        print("\n".join(list(data_index[0])), file=fout)
