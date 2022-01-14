import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import preprocess 
import pca_check 
from dlearn import ae
import importlib


df = pd.read_csv("../output/cd4_c13_c18_counts.txt.gz",sep="\t")
df = preprocess.filter_genes(df)

###
# pca_check.run_pca(df.iloc[:,1:df.shape[1]-1],df.iloc[:,-1].values)

###
device = ae.torch.device('cuda' if ae.torch.cuda.is_available() else 'cpu')
input_dims = df.iloc[:,1:-1].shape[1]
latent_dims = 32
autoencoder = ae.Autoencoder(input_dims,latent_dims).to(device)
y = pd.factorize(df.iloc[:,-1])[0]
data = ae.load_data(df.iloc[:,1:-1].to_numpy(),y,device)
autoencoder = ae.train(autoencoder,data,device)

# def plot_latent(autoencoder, data, num_batches=100):
# for i, (x,y) in enumerate(data):
# 	z = autoencoder.encoder(x.to(device))
# 	z = z.to('cpu').detach().numpy()
# 	plt.scatter(z[:, 0], z[:, 1],c=y, cmap='tab10')
# 	if i > 100:
# 		break
# plt.savefig("../output/ae.png");plt.close()


