import numpy as np
import pandas as pd
from collections import namedtuple
from gen_util.io import read_config
import matplotlib.pylab as plt
import umap
import umap.plot 

config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/scmetm.yaml"
params = read_config(config)
args = namedtuple('Struct',params.keys())(*params.values())

df_z = pd.read_csv(args.home+args.output+args.model["out"]+args.model["mfile"]+"etm_zz_data.tsv",sep="\t")

umap_2d = umap.UMAP(n_components=2, init='random', random_state=0)
proj_2d = umap_2d.fit(df_h.iloc[:,1:])
# plt.scatter(proj_2d[:,0],proj_2d[:,1],s=0.001)

df_h = pd.read_csv(args.home+args.output+args.model["out"]+args.model["mfile"]+"etm_hh_data.tsv",sep="\t")
df_h.columns = [ x.replace("hh","k") for x in df_h.columns]
topic_labels = df_h.iloc[:,1:].idxmax(axis=1)


from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=df_h.shape[1]-1, random_state=0).fit(df_h.iloc[:,1:].to_numpy())
umap.plot.points(proj_2d,labels=kmeans.labels_)

plt.savefig("umap_hh_v4.png",dpi=300);plt.close()

# ##test to get neighbors
# embeddings = umap.UMAP(n_neighbors=15, min_dist=0.5).fit_transform(df_h.iloc[:,1:])
# knn = umap.umap_.nearest_neighbors(embeddings, 
#         n_neighbors=15, metric='euclidean', 
#         metric_kwds={}, angular=True, 
#         random_state=np.random.RandomState(42))

import pickle

f_name = 'umap_model.pkl'
pickle.dump(umap_2d, open(f_name, 'wb'))