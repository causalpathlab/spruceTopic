import pandas as pd
import numpy as np
import preprocess 
import component_check 
from dlearn import ae
import importlib

data="cd4_cd8_500cells_per_tissue_counts.txt.gz"
df = pd.read_csv("../output/"+data,sep="\t")
df = df[ ["cell"]+\
        [x for x in df.columns if x not in["cell","sample"]]+\
        ["sample"]]
df = preprocess.filter_minimal(df)

y = np.array([1 if "CD4" in x else 0 for x in df.iloc[:,-1]])
# y = pd.factorize(df.iloc[:,-1])[0]

component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"tsne")

#### ae tests
device = ae.torch.device('cuda' if ae.torch.cuda.is_available() else 'cpu')
input_dims = df.iloc[:,1:-1].shape[1]
latent_dims = 32

autoencoder = ae.Autoencoder(input_dims,latent_dims).to(device)
data = ae.load_data(df.iloc[:,1:-1].to_numpy(),y,device)
autoencoder = ae.train(autoencoder,data,device)
ae.plot_ae_components(data,autoencoder,device,"cd4_cd8_ae")


df_z = pd.concat([df[["cell"]],ae.get_encoded_z(df,autoencoder,device)],axis=1)
df_z = df_z.set_index("cell")

from sklearn.neighbors import NearestNeighbors
import igraph

nbrs = NearestNeighbors(n_neighbors=10).fit(df_z)
nbrs_mat = nbrs.kneighbors_graph(df_z).toarray()
np.fill_diagonal(nbrs_mat, 0)

g = igraph.Graph.Adjacency((nbrs_mat > 0).tolist())
g.es['weight'] = nbrs_mat[nbrs_mat.nonzero()]
g.vs['label'] = df_z.index


visual_style = {}
# Define colors used for outdegree visualization
colours = ['#fecc5c', '#a31a1c']
# Set bbox and margin
visual_style["bbox"] = (3000,3000)
visual_style["margin"] = 17
# Set vertex colours
visual_style["vertex_color"] = 'grey'
# Set vertex size
visual_style["vertex_size"] = 20
# Set vertex lable size
visual_style["vertex_label_size"] = 0
# Don't curve the edges
visual_style["edge_curved"] = False
# Set the layout
my_layout = g.layout_fruchterman_reingold()
visual_style["layout"] = my_layout
# Plot the graph
igraph.plot(g, target='../output/interaction.pdf', **visual_style)
