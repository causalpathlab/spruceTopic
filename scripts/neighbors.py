import pandas as pd
import numpy as np
import random
import matplotlib.pylab as plt
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from itertools import combinations

def get_neighbors(df,nbr_mode):
    X = df.iloc[:,1:-1].values
    nbrs = NearestNeighbors(n_neighbors=nbr_mode[1]).fit(X)
    distances, indices = nbrs.kneighbors(X)
    cell_pairs = []

    if nbr_mode[0] == "one_cell":
        i = indices[0]
        cell_pair_indexes = [":".join(map(str, comb)) for comb in combinations(i, 2)]
        for cell_pair_index in cell_pair_indexes:
            i1,i2 = int(cell_pair_index.split(":")[0]),int(cell_pair_index.split(":")[1])
            cell_pair =  ":".join(df.cell[[i1,i2]].sort_values().values)
            if cell_pair not in cell_pairs:
                cell_pairs.append(cell_pair)
        return [x.split(':') for x in cell_pairs]

    elif nbr_mode[0] == "random_cell":
        for i in indices:
            cell_pair_indexes = [":".join(map(str, comb)) for comb in combinations(i, 2)]
            for cell_pair_index in cell_pair_indexes:
                i1,i2 = int(cell_pair_index.split(":")[0]),int(cell_pair_index.split(":")[1])
                cell_pair =  ":".join(df.cell[[i1,i2]].sort_values().values)
                if cell_pair not in cell_pairs:
                    cell_pairs.append(cell_pair)
        return [x.split(':') for x in cell_pairs]

def genes_from_ligand_receptor_db(mode):
    lr_db="../database/celltalkdb_v20220131_human_lr_pair.txt"
    df = pd.read_csv(lr_db,sep="\t")
    if mode == "gene":
        genes = []
        for gp in df["lr_pair"].values:
            gp = gp.split("_");genes.append(gp[0]);genes.append(gp[1])
        return genes
    elif mode == "lrpair":
        return df["lr_pair"].values

def keep_lr_genes():
    ### select only genes present in the ligand-receptor database 
    # need to run first time to filter out genes
    data = "cd4_cd8_500cells_per_tissue_counts.txt.gz"
    df_exp = pd.read_csv("../output/cell_data/"+data,sep="\t")
    lr_genes = genes_from_ligand_receptor_db("gene")
    selected_genes = [x for x in df_exp.columns if x in lr_genes ]
    df_exp = df_exp[["cell"]+selected_genes+["sample"]]
    df_exp.to_csv("cd4_cd8_500cells_per_tissue_counts_l_r_pair.txt.gz",index=False,sep="\t",compression="gzip")

def adhoc_check(df_exp,cell_pairs,genes):
    ave_exp= 5
    res ={}
    for indx,grp in enumerate(cell_pairs):
        cell_group_indexes = df_exp[df_exp["cell"].isin(grp)].index
        for gene_pair in genes:
            df_cell_group = df_exp.loc[cell_group_indexes,gene_pair.split("_")]
            if df_cell_group.mean()[0] >ave_exp and df_cell_group.mean()[1] >ave_exp:
                if gene_pair not in res:
                    res[gene_pair] = 1
                else:
                    res[gene_pair] += 1
                    print(gene_pair,res[gene_pair])
    df_res = pd.DataFrame([res])
    df_res = df_res.T
    df_res.to_csv("../output/cell_data/interaction_result.csv",index=False,sep="\t")

def interaction_method_gene(xij,m):
    if m == "add":
        return xij[0,:]+xij[1,:]
    elif m == "mult":
        return xij[0,:]*xij[1,:]
    elif m == "thr":
        xij_hat =[]
        cutoff=5
        for indx in range(xij.shape[1]) :
            if xij[0,indx]>cutoff and xij[1,indx]>cutoff:
                xij_hat.append(1)
            else:
                xij_hat.append(0)
        return xij_hat

def incidence_mat_gene(df_exp,cell_pairs,genes,interaction_mode):
    mat_xij = np.empty((0,len(genes)), int)
    for pair in cell_pairs:
        x_ij = df_exp.loc[df_exp["cell"].isin(pair),genes].values
        x_ij = interaction_method_gene(x_ij,interaction_mode)
        mat_xij = np.append(mat_xij, np.array([x_ij]), axis=0)
    return np.matrix(mat_xij)

def interaction_method_lrpair(x,m):
    if m == "ep": #expression product
        return np.max( [ (x[0,0] * x[1,1]) - (x[0,1] * x[1,0]) , \
                         (x[0,1] * x[1,0]) - (x[0,0] - x[1,1]) ] )
    elif m== "et": # expression thresholding
        thr = 3
        if (x[0,0] > 1 and x[1,1] > 1) or (x[1,0] > 1 and x[0,1] > 1):
            return 1
        else:
            return 0

def incidence_mat_lrpair(df_exp,cell_pairs,lr_pairs):
    mat_xij = np.empty((0,len(lr_pairs)), int)
    for cell_pair in cell_pairs:
        x_ij = []
        for lr_pair in lr_pairs:
            cell_lr_table = df_exp.loc[df_exp["cell"].isin(cell_pair),lr_pair.split('_')].values
            x_ij.append(interaction_method_lrpair(cell_lr_table,"ep"))
        mat_xij = np.append(mat_xij, np.array([x_ij]), axis=0)
    return np.matrix(mat_xij)

