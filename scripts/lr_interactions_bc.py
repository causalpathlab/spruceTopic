import pandas as pd
import neighbors as nbr
import lr_map as lrm
import preprocess
import importlib
import random
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
plt.rcParams["figure.figsize"] = [15.50, 3.50]
plt.rcParams["figure.autolayout"] = True

### get lr genes
genes = lrm.genes_from_lr_db("gene")

###bc data
exp_data = "bc_counts.txt.gz"
df_exp = pd.read_csv("../output/bc_sc_data/"+exp_data,sep="\t")
keep_lr_genes = ["cell"]+[x for x in df_exp.columns if x in genes ] +["sample"]
df_exp = df_exp[keep_lr_genes]
### get z from etm 
zz="../output/bc_sc_data/sc_zz_sc_epochs_2000_data.csv"
df = pd.read_csv(zz,sep="\t")
df["cell"] =df_exp["cell"]
df["sample"] =df_exp["sample"]
df = df[["cell"]+[x for x in df.columns[1:-2] ] +["sample"]]


##normalize
df_exp = df_exp.loc[:, (df_exp != 0).any(axis=0)]
row_sum = df_exp.iloc[:,1:-1].sum(1).values
df_exp.iloc[:,1:-1] = df_exp.iloc[:,1:-1].apply(lambda x: x/row_sum, axis=0)

# get genes from lr db - remove genes present in lr db but not in experiment data
genes = [ x for x in genes if x in df_exp.columns]
genes = pd.Series(genes).unique()

### get graph 
g = lrm.graph_lr_db(genes)

####
# nbr_mode = ("oncp",15)
nbr_mode = ("ancp",2)
# nbr_mode = ("rncp",2)

#### interactions
cell_pairs_batch = nbr.get_neighbors(df,nbr_mode)
mat_xij = nbr.incidence_mat_lrgraph(g,cell_pairs_batch,df_exp)

df_meta = pd.read_csv("../input/bc_data/GSE75688_final_sample_information.txt",sep="\t")
meta = {}
for x,y in zip(df_meta["sample"],df_meta["index2"]):meta[x]=y 

cell_pairs_type = []
for x in cell_pairs_batch:
    try:
        id = x[0]+'_'+ meta[x[0]]+'/'+ x[1]+'_'+ meta[x[1]]
    except:
        id = x[0]+'/'+ x[1]
    cell_pairs_type.append(id)


df_mat = pd.DataFrame(mat_xij)
df_mat.columns = [g.vs[e.source]['name']+'_'+g.vs[e.target]['name'] for e in g.es]
df_mat.index = cell_pairs_type

tumor_immune = []
for i in cell_pairs_type:
    if "Tumor" in i and "Stromal" in i:tumor_immune.append(i)
df_mat = df_mat[df_mat.index.isin(tumor_immune)]

plt.plot(df_mat.values.flatten())
plt.savefig("../df_vals.png");plt.close()

noise_thr = 0.00005
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

