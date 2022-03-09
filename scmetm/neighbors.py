import pandas as pd
import numpy as np
import random
from sklearn.neighbors import NearestNeighbors
from itertools import combinations

def get_neighbors(df,nbr_mode):
    X = df.iloc[:,1:-1].values
    nbrs = NearestNeighbors(n_neighbors=nbr_mode[1]).fit(X)
    distances, indices = nbrs.kneighbors(X)
    cell_pairs = []

    if nbr_mode[0] == "oncp":
        i = indices[np.random.randint(len(indices), size=1)][0]
        # i = indices[0]
        cell_pair_indexes = [comb for comb in combinations(i, 2)]
        for cell_pair_index in cell_pair_indexes:
            i1,i2 = int(cell_pair_index[0]),int(cell_pair_index[1])
            cell_pair =  list(df.cell[[i1,i2]].sort_values().values)
            if cell_pair not in cell_pairs:
                cell_pairs.append(cell_pair)
        return cell_pairs

    elif nbr_mode[0] == "ancp":
        for i in indices:
            cell_pair_indexes = [ comb for comb in combinations(i, 2)]
            for cell_pair_index in cell_pair_indexes:
                i1,i2 = int(cell_pair_index[0]),int(cell_pair_index[1])
                cell_pair =  list(df.cell[[i1,i2]].sort_values().values)
                if cell_pair not in cell_pairs:
                    cell_pairs.append(cell_pair)
        return cell_pairs

    elif nbr_mode[0] == "rncp":
        selected_pair = []
        discarded_cells = []
        for indx in indices:
            i,j = indx[0],indx[1]
            if i not in discarded_cells and j not in discarded_cells:
                selected_pair.append((i,j))
            for x in indx:discarded_cells.append(x)
        for cell_pair_index in selected_pair:
            i1,i2 = int(cell_pair_index[0]),int(cell_pair_index[1])
            cell_pair =  list(df.cell[[i1,i2]].sort_values().values)
            cell_pairs.append(cell_pair)
        return cell_pairs

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
    return df_res
    
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
        # return np.max( [ (x[0,0] * x[1,1]) , (x[0,1] * x[1,0]) ] )
        return np.max( [ (x[0,0] * x[1,1]) - (x[0,1] * x[1,0]) , \
                         (x[0,1] * x[1,0]) - (x[0,0] * x[1,1]) ] )
    elif m== "et": # expression thresholding
        thr = 3
        if (x[0,0] > thr and x[1,1] > thr) or (x[1,0] > thr and x[0,1] > thr):
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

def incidence_mat_lrgraph(g,cell_pairs,df_exp):
    mat_xij = np.empty((0,len(g.es)), int)
    for cell_pair in cell_pairs:
        x_ij = []
        for e in g.es:
            g1 = g.vs[e.source]['name']
            g2 = g.vs[e.target]['name']
            cell_lr_table = df_exp.loc[df_exp["cell"].isin(cell_pair),[g1,g2]].values
            x_ij.append(interaction_method_lrpair(cell_lr_table,"ep"))
        mat_xij = np.append(mat_xij, np.array([x_ij]), axis=0)
    return np.matrix(mat_xij)
