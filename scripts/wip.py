import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import preprocess 
import component_check 
from dlearn import ae
import importlib

# data="cd4_c13_c18_counts.txt.gz"
# df = pd.read_csv("../output/"+data,sep="\t")
# df = preprocess.filter_genes(df)
# df = df[(df.cluster=="BC_CD4_c13")|(df.cluster=="BC_CD4_c18")]

data="pbmc3k_scanpy_filtered_counts.txt.gz"
df = pd.read_csv("../output/"+data)
df = df.reset_index()
df["sample"]="c13"

y = np.array([1 if "c13" in x else 0 for x in df.iloc[:,-1]])
component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"tsne")

#### ae tests
device = ae.torch.device('cuda' if ae.torch.cuda.is_available() else 'cpu')
input_dims = df.iloc[:,1:-1].shape[1]
latent_dims = 32

autoencoder = ae.Autoencoder(input_dims,latent_dims).to(device)
data = ae.load_data(df.iloc[:,1:-1].to_numpy(),y,device)
autoencoder = ae.train(autoencoder,data,device)

autoencoder = ae.VariationalAutoencoder(input_dims,latent_dims).to(device)
data = ae.load_data(df.iloc[:,1:-1].to_numpy(),y,device)
autoencoder = ae.vaetrain(autoencoder,data,device)


for i, (x,y) in enumerate(data):
    z = autoencoder.encoder(x.to(device))
    z = z.to('cpu').detach().numpy()
    plt.scatter(z[:, 0], z[:, 1],c=y, cmap='tab10')
plt.savefig("../output/pbmc_ae.png");plt.close()



