import pandas as pd
import re
import sys
import os

blast_file = sys.argv[1]
sample=sys.argv[2]
data = pd.read_table(blast_file, header=None)
all_se1_ss2_se2_ss3 = {(data[9][i], data[13][i], data[14][i], data[18][i]) for i in range(len(data))}
if not os.path.exists("reads_group"):
    os.mkdir("reads_group")
for item in list(all_se1_ss2_se2_ss3):
    data_tmp = data[(data[9] == item[0]) & (data[13] == item[1]) & (data[14] == item[2]) & (data[18] == item[3])]
    data_tmp.to_csv("reads_group/" + str(item[0]) + "_" + str(item[1]) + "_" + str(item[2]) + "_" + str(item[3]) + "_" + blast_file, header=False, index=False, sep="\t")
    repsize_1 = list(set(data_tmp[3]))
    repsize_1.sort()
    repsize_str_list_1 = [str(item) for item in repsize_1]
    repsize_str_1 = ",".join(repsize_str_list_1)
    repsize_2 = list(set(data_tmp[4]))
    repsize_2.sort()
    repsize_str_list_2 = [str(item) for item in repsize_2]
    repsize_str_2 = ",".join(repsize_str_list_2)
    with open("count_se1_ss2_se2_ss3_" + blast_file, "at") as f_count:
        print("\t".join([str(item[0]), str(item[1]), str(item[2]), str(item[3]), str(len(data_tmp)), repsize_str_1, repsize_str_2, sample]), file=f_count)

data_count = pd.read_table("count_se1_ss2_se2_ss3_" + blast_file, header=None)
all_se1_ss2_se2_ss3 = {(data_count[0][i], data_count[1][i], data_count[2][i], data_count[3][i]) for i in range(len(data_count))}
all_list = list(all_se1_ss2_se2_ss3)
all_list.sort(key=lambda x:x[0])
df_count=pd.DataFrame(columns=["se1", "ss2", "se2", "ss3", sample], index=all_list)
for i in range(len(data_count)):
    index_tuple = (data_count[0][i], data_count[1][i], data_count[2][i], data_count[3][i])
    df_count[data_count[7][i]][index_tuple] = data_count[4][i]
df_count.fillna(0, inplace=True)
df_count.to_excel(sample + "_se1_ss2_se2_ss3_annotation.xlsx")
df_count = pd.read_excel(sample + "_se1_ss2_se2_ss3_annotation.xlsx")
df_count.rename(columns={"Unnamed: 0":"se1_ss2_se2_ss3"}, inplace=True)
data1 = pd.read_table("count_se1_ss2_se2_ss3_" + blast_file, header=None)
for i in range(len(df_count)):
    cords = re.findall("\d+", df_count.loc[i, "se1_ss2_se2_ss3"])
    data1_tmp = data1[(data1[0] == int(cords[0])) & (data1[1] == int(cords[1])) & (data1[2] == int(cords[2])) & (data1[3] == int(cords[3]))]
    df_count.loc[i, "se1"] = int(cords[0])
    df_count.loc[i, "ss2"] = int(cords[1])
    df_count.loc[i, "ss2"] = int(cords[2])
    df_count.loc[i, "se3"] = int(cords[3])    
    df_count.loc[i, "index"] = i + 1
    df_count.loc[i, "repsize_1"] = str(data1_tmp.iloc[0].iat[5]).lstrip()
    df_count.loc[i, "min_repsize_1"] = int(re.findall("\d+", df_count.loc[i, "repsize_1"])[0])
    df_count.loc[i, "max_repsize_1"] = int(re.findall("\d+", df_count.loc[i, "repsize_1"])[-1])
    df_count.loc[i, "mid_repsize_1"] = (df_count.loc[i, "min_repsize_1"] + df_count.loc[i, "max_repsize_1"])/2
    df_count.loc[i, "repsize_2"] = str(data1_tmp.iloc[0].iat[6]).lstrip()
    df_count.loc[i, "min_repsize_2"] = int(re.findall("\d+", df_count.loc[i, "repsize_2"])[0])
    df_count.loc[i, "max_repsize_2"] = int(re.findall("\d+", df_count.loc[i, "repsize_2"])[-1])
    df_count.loc[i, "mid_repsize_2"] = (df_count.loc[i, "min_repsize_2"] + df_count.loc[i, "max_repsize_2"])/2

df_count = df_count[["index", "se1_ss2_se2_ss3", "se1", "ss2", "se2", "ss3", "repsize_1", "min_repsize_1", "mid_repsize_1", "max_repsize_1", "repsize_2", "min_repsize_2", "mid_repsize_2", "max_repsize_2", sample]]
df_count.sort_values("mid_repsize_1", inplace = True, ascending=False)
df_count.to_excel(sample + "_se1_ss2_se2_ss3_annotation.xlsx", index=False)
df_count.to_excel(sample + "_se1_ss2_se2_ss3_backup.xlsx", index=False)
if not os.path.exists("reads_group_ids"):
    os.mkdir("reads_group_ids")
for i in range(len(df_count)):
    data_index = data[(data[9] == item[0]) & (data[13] == item[1]) & (data[14] == item[2]) & (data[18] == item[3])]
    with open("reads_group_ids/" + str(int(df_count.loc[i, "index"])) + "_ids.list", "wt") as fout:
        print("\n".join(list(data_index[0])), file=fout)
