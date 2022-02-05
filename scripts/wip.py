import pandas as pd
import neighbors as nbr
import importlib
import random
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
plt.rcParams["figure.figsize"] = [15.50, 3.50]
plt.rcParams["figure.autolayout"] = True

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
cell_pairs = nbr.get_neighbors(df)

N=100
cell_pairs_batch = random.sample(cell_pairs,N)

cell_pairs = cell_pairs_batch

mode="thr"
mat_xij = nbr.incidence_mat_gene(df_exp,cell_pairs_batch,genes,mode)

cell_pairs_type = []
for cp in cell_pairs_batch:
    types = df_exp.loc[df_exp["cell"].isin(cp),["sample"]].values;cell_pairs_type.append(types.flatten()[0]+"/"+types.flatten()[1])


df_mat = pd.DataFrame(mat_xij)
df_mat.columns = genes

plt.plot(df_mat.values.flatten())
plt.savefig("../df_vals.png");plt.close()


## check k-means
from sklearn.cluster import KMeans
df_mat = df_mat.T
kmeans = KMeans(n_clusters=5, random_state=0).fit(df_mat.values)
df_mat["cluster"] = kmeans.labels_
df_mat = df_mat.sort_values("cluster")
df_mat["cluster"].value_counts()

df_mat = df_mat[df_mat.cluster!=0]


sns.heatmap(df_mat.T.iloc[0:-1,:],cmap="binary")
plt.savefig("../test.png");plt.close()

## need to figure out how to plot edges with cluster color
import lr_map 
import igraph

colors = [plt.cm.Pastel1(x) for x in range(5)]
color_dict = dict(zip(range(5), colors))
edge_colors = [color_dict[x] for x in kmeans.labels_]


g = lr_map.plot_df(df_mat)


g.vs["name"] = cell_pairs_type


to_delete_ids = [v.index for v in g.vs if v.degree()<1]
g.delete_vertices(to_delete_ids)
g.es['color'] = edge_colors
visual_style={}
visual_style["vertex_label"] = g.vs["name"]
visual_style["vertex_size"] = 10
visual_style["vertex_label_size"] = 8
visual_style["vertex_label_color"] = "darkblue"
visual_style["edge_width"] = g.es['weight']
visual_style["layout"] = g.layout_fruchterman_reingold()
igraph.plot(g, target="../network.pdf",**visual_style)

