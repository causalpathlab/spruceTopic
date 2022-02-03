import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from itertools import combinations

def get_neighbors(df):
    X = df.iloc[:,1:-1].values
    nbrs = NearestNeighbors(n_neighbors=5).fit(X)
    distances, indices = nbrs.kneighbors(X)
    cell_pairs = []
    for i in indices:
        cell_pair_indexes = [":".join(map(str, comb)) for comb in combinations(i, 2)]
        for cell_pair_index in cell_pair_indexes:
            i1,i2 = int(cell_pair_index.split(":")[0]),int(cell_pair_index.split(":")[1])
            cell_pair =  ":".join(df.cell[[i1,i2]].sort_values().values)
            if cell_pair not in cell_pairs:
                cell_pairs.append(cell_pair)
    return [x.split(':') for x in cell_pairs]

def genes_from_ligand_receptor_db():
    lr_db="../database/celltalkdb_v20220131_human_lr_pair.txt"
    df = pd.read_csv(lr_db,sep="\t")
    genes = []
    for gp in df["lr_pair"].values:
        gp = gp.split("_");genes.append(gp[0]);genes.append(gp[1])
    return genes

def keep_lr_genes():
    ### select only genes present in the ligand-receptor database 
    # need to run first time to filter out genes
    data = "cd4_cd8_500cells_per_tissue_counts.txt.gz"
    df_exp = pd.read_csv("../output/cell_data/"+data,sep="\t")
    lr_genes = genes_from_ligand_receptor_db()
    selected_genes = [x for x in df_exp.columns if x in lr_genes ]
    df_exp = df_exp[["cell"]+selected_genes+["sample"]]
    df_exp.to_csv("cd4_cd8_500cells_per_tissue_counts.txt.gz",index=False,sep="\t",compression="gzip")

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

def interaction_method(xij,m):
    if m == "add":
        return xij[0,:]+xij[1,:]
    elif m == "mult":
        return xij[0,:]*xij[1,:]
    elif m == "thr":
        xij_hat =[]
        cutoff=1
        for indx in range(xij.shape[1]) :
            if xij[0,indx]>cutoff & xij[1,indx]>cutoff:
                xij_hat.append(1)
            else:
                xij_hat.append(0)
        return xij_hat

def incidence_mat(df_exp,cell_pairs,genes,interaction_mode):
    mat_xij = np.empty((0,len(genes)), int)
    for pair in cell_pairs:
        x_ij = df_exp.loc[df_exp["cell"].isin(pair),genes].values
        x_ij = interaction_method(x_ij,interaction_mode)
        mat_xij = np.append(mat_xij, np.array([x_ij]), axis=0)
    return np.matrix(mat_xij)

######
data = "cd4_cd8_500cells_per_tissue_counts_l_r_pair.txt.gz"
df_exp = pd.read_csv("../output/cell_data/"+data,sep="\t")

# get genes from lr db
# remove genes present in lr db but not in experiment data
genes = genes_from_ligand_receptor_db()
check_genes = [x for x in genes if x not in df_exp.columns]
genes = [ x for x in genes if x not in check_genes]
genes = pd.Series(genes).unique()

### get z from etm 
zz="../output/cell_data/sc_zz_sc_epochs_2000_data.csv"
df = pd.read_csv(zz,sep="\t")
cell_pairs = get_neighbors(df)

interaction_mode="add"
mat_xij = incidence_mat(df_exp,cell_pairs,genes,interaction_mode)
df_mat = pd.DataFrame(mat_xij)

sns.set()
sns.heatmap(df_mat,cmap="YlGnBu")
plt.savefig("../test.png");plt.close()