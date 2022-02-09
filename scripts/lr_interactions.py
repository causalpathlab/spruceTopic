import pandas as pd
import neighbors as nbr
import lr_map as lrm
import importlib
import random
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
plt.rcParams["figure.figsize"] = [15.50, 3.50]
plt.rcParams["figure.autolayout"] = True

######
exp_data = "cd4_cd8_500cells_per_tissue_counts_l_r_pair.txt.gz"
df_exp = pd.read_csv("../output/cell_data/"+exp_data,sep="\t")
##normalize
df_exp = df_exp.loc[:, (df_exp != 0).any(axis=0)]
row_sum = df_exp.iloc[:,1:-1].sum(1).values
df_exp.iloc[:,1:-1] = df_exp.iloc[:,1:-1].apply(lambda x: x/row_sum, axis=0)


# get genes from lr db - remove genes present in lr db but not in experiment data
genes = lrm.genes_from_lr_db("gene")
genes = [ x for x in genes if x in df_exp.columns]
genes = pd.Series(genes).unique()


### get graph 
g = lrm.graph_lr_db(genes)

### get z from etm 
zz="../output/cell_data/sc_zz_sc_epochs_2000_data.csv"
df = pd.read_csv(zz,sep="\t")

####
'''
neighbor cell search mode:
    oncp - one neighbor cell pair, select any one cell and get its neighboring cells 
    ancp - all neighbor cell pair, get neighboring cell for all cells
    rncp - random neighbor cell pair, get non-overlaping neighboring cell pair
'''
nbr_mode = ("oncp",25)
# nbr_mode = ("rncp",25)

#### interactions
cell_pairs = nbr.get_neighbors(df,nbr_mode)
N=25
cell_pairs_batch = random.sample(cell_pairs,N)
mat_xij = nbr.incidence_mat_lrgraph(g,cell_pairs_batch,df_exp)
cell_pairs_type = []
for cp in cell_pairs_batch:
    types = df_exp.loc[df_exp["cell"].isin(cp),["sample"]].values;cell_pairs_type.append(types.flatten()[0]+"/"+types.flatten()[1])

df_mat = pd.DataFrame(mat_xij)
df_mat.columns = [g.vs[e.source]['name']+'_'+g.vs[e.target]['name'] for e in g.es]
df_mat.index = cell_pairs_type

# plt.plot(df_mat.values.flatten())
# plt.savefig("../df_vals.png");plt.close()

noise_thr = 0.01
df_mat = df_mat.loc[:, (df_mat > noise_thr).any(axis=0)]

## cluster edges k-means
from sklearn.cluster import KMeans
df_mat = df_mat.T
kmeans = KMeans(n_clusters=5, random_state=0).fit(df_mat.values)
df_mat["cluster"] = kmeans.labels_
df_mat = df_mat.sort_values("cluster")
df_mat["cluster"].value_counts()
plt.rcParams["figure.figsize"] = [15.50, 10.50]
plt.rcParams["figure.autolayout"] = True
sns.heatmap(df_mat.T.iloc[0:-1,:],cmap="YlGnBu")
plt.savefig("../output/nbr_cell_lr_db_heatmap_"+nbr_mode[0]+"_"+str(nbr_mode[1])+".png");plt.close()

g = lrm.update_graph(g,df_mat)
lrm.plot_graph(g,df_mat,nbr_mode[0]+"_"+str(nbr_mode[1]))

