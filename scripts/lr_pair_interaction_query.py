import pandas as pd
import neighbors as nbr
import importlib
import random
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np

######
data = "cd4_cd8_500cells_per_tissue_counts_l_r_pair.txt.gz"
df_exp = pd.read_csv("../output/cell_data/"+data,sep="\t")

# get genes from lr db
# remove genes present in lr db but not in experiment data
genes = nbr.genes_from_ligand_receptor_db("gene")
lr_pairs = nbr.genes_from_ligand_receptor_db("lrpair")
check_genes = [x for x in genes if x not in df_exp.columns]
genes = [ x for x in genes if x not in check_genes]
genes = pd.Series(genes).unique()
lr_pairs = [ x for x in lr_pairs if x.split('_')[0] not in check_genes and x.split('_')[1] not in check_genes  ]

### get z from etm 
zz="../output/cell_data/sc_zz_sc_epochs_2000_data.csv"
df = pd.read_csv(zz,sep="\t")

####
# nbr_mode = ("one_cell",25)
nbr_mode = ("random_cell",2)


####
cell_pairs = nbr.get_neighbors(df,nbr_mode)
print("number of cell pair selected:")
print(len(cell_pairs))

##
N=25
cell_pairs_batch = random.sample(cell_pairs,N)
mat_xij = nbr.incidence_mat_lrpair(df_exp,cell_pairs_batch,lr_pairs)

cell_pairs_type = []
for indx,cp in enumerate(cell_pairs_batch):
    types = df_exp.loc[df_exp["cell"].isin(cp),["sample"]].values;cell_pairs_type.append(str(indx)+"_"+types.flatten()[0]+"/"+types.flatten()[1])

df_mat = pd.DataFrame(mat_xij)
df_mat.columns = lr_pairs
df_mat.index = cell_pairs_type

noise_thr = 10
df_mat = df_mat.loc[:, (df_mat > noise_thr).any(axis=0)]


## check k-means
from sklearn.cluster import KMeans
df_mat = df_mat.T
kmeans = KMeans(n_clusters=5, random_state=0).fit(df_mat.values)
df_mat["cluster"] = kmeans.labels_
df_mat = df_mat.sort_values("cluster")
df_mat["cluster"].value_counts()
plt.rcParams["figure.figsize"] = [15.50, 10.50]
plt.rcParams["figure.autolayout"] = True
sns.heatmap(df_mat.T.iloc[0:-1,:],cmap="YlGnBu")
plt.savefig("../output/nbr_cell_lr_pairs_heatmap_"+nbr_mode[0]+"_"+str(nbr_mode[1])+".png");plt.close()


