# =-*- coding: utf-8 -*-
import argparse
import numpy as np
import pandas as pd
import numba as nb

# Input parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Process command-line arguments.")
    parser.add_argument("-s", type=str, required=True, help="StringDB dataframe input")
    parser.add_argument("-sl", type=int, default=450, help="StringDB score limit")
    parser.add_argument(
        "-dfexpr", type=str, required=True, help="Expression dataframe input"
    )
    parser.add_argument(
        "-dfanno", type=str, required=True, help="Annotation dataframe input"
    )
    parser.add_argument(
        "-dic", type=str, required=True, help="Dictionary for ENSP to symbol"
    )
    parser.add_argument("-o", type=str, default="", help="Output name")
    return parser.parse_args()

args = parse_arguments()

# Initialize variables
stringdb = pd.read_csv(args.s, delimiter=" ")
string_limit = min(max(args.sl, 0), 999)
expr = pd.read_csv(args.dfexpr)
anno = pd.read_csv(args.dfanno)
dict = pd.read_csv(args.dic)
output_name = str(args.o)

# Stringdb preprocessing
stringdb["protein1"] = stringdb["protein1"].str.slice(start=5)
stringdb["protein2"] = stringdb["protein2"].str.slice(start=5)
stringdb = stringdb[stringdb["combined_score"] >= string_limit]
stringdb = stringdb.reset_index(drop=True)

# Construct dictionary
length = len(np.unique(dict["symbol"].values))

d = {dict["ensp"][i]: int(dict["id"][i]) for i in range(len(dict))}
d2 = {int(dict["id"][i]): dict["symbol"][i] for i in range(len(dict))}
string_bool = np.zeros([length, length])

for i in range(len(stringdb)):
    if (stringdb["protein1"][i] in d) & (stringdb["protein2"][i] in d):
        string_bool[d[stringdb["protein1"][i]], d[stringdb["protein2"][i]]] = 1

np.fill_diagonal(string_bool, 0)
string_bool = string_bool == 1

# Reference sample and purturbed sample
purturbed = anno[anno["time"] != "Baseline"]["geo"]
reference = anno[anno["time"] == "Baseline"]["geo"]

# Entropy calculation
@nb.njit(parallel=False)
def entropy(pc, sd, idlist1, idlist2, edge_entropy, edge_sd):
    for i in range(len(idlist1)):
        p1 = idlist1[i]
        p2 = idlist2[i]

        left_not_zero = np.where(pc[p1] != 0)[0]
        right_not_zero = np.where(pc[p2] != 0)[0]

        left_not_zero = np.delete(left_not_zero, np.where(left_not_zero == p2)[0])
        right_not_zero = np.delete(right_not_zero, np.where(right_not_zero == p1)[0])

        len_left_not_zero = len(left_not_zero)
        len_right_not_zero = len(right_not_zero)

        if len_left_not_zero < 2 and len_right_not_zero < 2:
            edge_entropy[i] = 0
            edge_sd[i] = 0
            continue

        left_weight = len_left_not_zero if len_left_not_zero > 0 else 0
        right_weight = len_right_not_zero if len_right_not_zero > 0 else 0

        if left_weight > 0 and right_weight > 0:
            left_prob = (pc[p1][left_not_zero]) / np.sum(pc[p1][left_not_zero])
            right_prob = (pc[p2][right_not_zero]) / np.sum(pc[p2][right_not_zero])

            left_entropy = -np.sum(left_prob * np.log2(left_prob))
            right_entropy = -np.sum(right_prob * np.log2(right_prob))

            edge_entropy[i] = (left_entropy * left_weight + right_entropy * right_weight) / (left_weight + right_weight)
            edge_sd[i] = (sd[p1] * left_weight + sd[p2] * right_weight) / (left_weight + right_weight)
        else:
            edge_entropy[i] = 0
            edge_sd[i] = 0
            

# Reference samples calculation
origin_frame = expr[reference].to_numpy()
origin_sd = np.std(origin_frame, ddof=1, axis=1)
origin_pc = np.corrcoef(origin_frame)
origin_pc = np.abs(origin_pc)
origin_pc = origin_pc * string_bool
edge_origin_list = []
num = len(origin_frame)
idlist1 = []
idlist2 = []
for i in range(num):
    for j in range(i - 1):
        if string_bool[i, j] != 0:
            edge_origin_list.append([d2[i], d2[j]])
            idlist1.append(i)
            idlist2.append(j)

edge_origin_entropy = np.zeros((len(idlist1)), dtype=np.float64)
edge_origin_sd = np.zeros((len(idlist1)), dtype=np.float64)
idlist1 = nb.typed.List(idlist1)
idlist2 = nb.typed.List(idlist2)

entropy(origin_pc, origin_sd, idlist1, idlist2, edge_origin_entropy, edge_origin_sd)

# Purturbed samples calculation
landscape = pd.DataFrame()
single_sample = expr[purturbed]
num_sample = single_sample.columns
num_gene = len(single_sample.iloc[:, 0])

for k, item in enumerate(num_sample):
    print(k, end="\n")
    new = np.hstack(
        (origin_frame, np.reshape(single_sample.iloc[:, k].values, (num_gene, 1)))
    )

    new_pc = np.corrcoef(new)
    new_pc = np.abs(new_pc)
    new_pc = new_pc * string_bool
    new_sd = np.std(new, ddof=1, axis=1)

    edge_append_entropy = np.zeros((len(idlist1)))
    edge_append_sd = np.zeros((len(idlist1)))

    entropy(new_pc, new_sd, idlist1, idlist2, edge_append_entropy, edge_append_sd)

    edge_append_entropy = np.array(edge_append_entropy)
    edge_append_entropy = np.abs(edge_append_entropy - edge_origin_entropy)
    edge_append_sd = np.array(edge_append_sd)
    edge_append_sd = np.abs(edge_append_sd - edge_origin_sd)

    landscape_pros = pd.DataFrame(edge_append_sd * edge_append_entropy)
    landscape_pros.columns = [num_sample[k]]
    landscape = pd.concat([landscape, landscape_pros], axis=1)

# Output
landscape = pd.concat(
    [pd.DataFrame(edge_origin_list, columns=["node1", "node2"]), landscape], axis=1
)
landscape = landscape.fillna(0)
landscape.to_csv("NIEE_" + str(output_name) + str(string_limit) + ".csv", index=False)
