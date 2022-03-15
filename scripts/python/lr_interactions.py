import sys
import pandas as pd
from scmetm import neighbors as nbr ,lr_map as lrm, preprocess
import random
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
from gen_util.io import read_config
from collections import namedtuple
import logging


config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/scmetm.yaml"
params = read_config(config)
args = namedtuple('Struct',params.keys())(*params.values())

logging.basicConfig(filename=args.home+args.output+args.model["out"]+args.model["mfile"]+".log",
						format='%(asctime)s %(levelname)-8s %(message)s',
						level=logging.INFO,
						datefmt='%Y-%m-%d %H:%M:%S')


# ######
df_exp = pd.read_pickle(args.home+args.input+args.raw_lr_data)
# df = pd.read_csv(zz,sep="\t")
df_h = pd.read_csv(args.home+args.output+args.model["out"]+args.model["mfile"]+"etm_hh_data.csv",sep="\t")


### get lr genes
genes = lrm.genes_from_lr_db("gene",args)
# # get genes from lr db - remove genes present in lr db but not in experiment data
genes = pd.Series(genes).unique()
genes = [ x for x in genes if x in df_exp.columns]

# ####
# '''
# neighbor cell search mode:
#     oncp - one neighbor cell pair, select any one cell and get its neighboring cells 
#     ancp - all neighbor cell pair, get neighboring cell for all cells
#     rncp - random neighbor cell pair, get non-overlaping neighboring cell pair
# '''
# nbr_mode = ("oncp",15)
nbr_mode = ("rncp",15)

# #### interactions
cell_pairs = nbr.get_neighbors(df_h,nbr_mode)
# N=25
# cell_pairs_batch = random.sample(cell_pairs,N)

# ### get graph 
g = lrm.graph_lr_db(genes,args)
mat_xij = nbr.incidence_mat_lrgraph(g,cell_pairs,df_exp)

# ##get cell meta data -
# cell_pairs_type = preprocess.cellid_to_meta(cell_pairs)


df_mat = pd.DataFrame(mat_xij)
df_mat.columns = [g.vs[e.source]['name']+'_'+g.vs[e.target]['name'] for e in g.es]
# df_mat.index = cell_pairs_type


plt.plot(df_mat.values.flatten())
plt.savefig("df_vals.png");plt.close()

sns.displot( x=df_mat.values.flatten(), kind="kde")
plt.savefig("df_vals_dist.png");plt.close()

noise_thr = 100
df_mat = df_mat.loc[:, (df_mat > noise_thr).any(axis=0)]

## cluster edges k-means
from sklearn.cluster import KMeans
df_mat = df_mat.T
kmeans = KMeans(n_clusters=3, random_state=0).fit(df_mat.values)
df_mat["cluster"] = kmeans.labels_
df_mat = df_mat.sort_values("cluster")
df_mat["cluster"].value_counts()
plt.rcParams["figure.figsize"] = [15.50, 10.50]
plt.rcParams["figure.autolayout"] = True
sns.heatmap(df_mat.T.iloc[0:-1,:],cmap="YlGnBu")
plt.savefig("nbr_cell_lr_db_heatmap_"+nbr_mode[0]+"_"+str(nbr_mode[1])+".png");plt.close()

g = lrm.update_graph(g,df_mat)
lrm.plot_graph(g,df_mat,nbr_mode[0]+"_"+str(nbr_mode[1]))

