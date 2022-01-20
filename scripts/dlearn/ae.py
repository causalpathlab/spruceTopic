
import torch; torch.manual_seed(0)
import torchvision
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

class Stacklayers(nn.Module):
	def __init__(self,input_size,layers_design):
		super(Stacklayers, self).__init__()

		self.layers = nn.ModuleList()
		self.input_size = input_size
		for next_size,activfunc in layers_design:
			self.layers.append(nn.Linear(input_size,next_size))
			self.layers.append(nn.BatchNorm1d(next_size))
			input_size = next_size
			self.layers.append(self.get_activation(activfunc))
		
	def forward(self, input_data):
		for layer in self.layers:
			input_data = layer(input_data)
		return input_data
	
	def get_activation(self, act):
		if act == 'tanh':
			act = nn.Tanh()
		elif act == 'relu':
			act = nn.ReLU()
		else:
			print('Defaulting to tanh activations...')
			act = nn.Tanh()
		return act 

class Decoder(nn.Module):
	def __init__(self, input_dims, latent_dims):
		super(Decoder, self).__init__()
		self.linear1 = nn.Linear(latent_dims, 128)
		self.linear2 = nn.Linear(128, input_dims)
	def forward(self, z):
		output = F.relu(self.linear1(z))
		x_hat = F.relu(self.linear2(output))
		return x_hat


class Autoencoder(nn.Module):
	def __init__(self,input_dims,latent_dims):
		super(Autoencoder, self).__init__()
		self.encoder = Stacklayers(input_dims, \
						[(128,"relu"),\
						(latent_dims,"relu")])
		self.decoder = Decoder(input_dims,latent_dims)
	def forward(self, x):
		z = self.encoder(x)
		return self.decoder(z)


class TabularDataset(Dataset):
	def __init__(self, x,y):
		self.X = x
		self.y = y	
	def __len__(self):
		return len(self.X)
	def __getitem__(self, idx):
		return [self.X[idx], self.y[idx]]

def load_data(x,y,device):
	x = x.astype('float32')
	x = torch.from_numpy(x).to(device)
	return DataLoader(TabularDataset(x,y), batch_size=64, shuffle=True)

def train(autoencoder, data,device, epochs=300):
	opt = torch.optim.Adam(autoencoder.parameters())
	mse = nn.MSELoss()
	for epoch in range(epochs):
		loss = 0
		for indx,(x,y) in enumerate(data):
			x = x.to(device)
			opt.zero_grad()
			x_hat = autoencoder(x)
			train_loss = mse(x_hat,x)
			train_loss.backward()
			opt.step()
			loss += train_loss.item()
		if epoch % 30 == 0:  
			print('====> Epoch: {} Average loss: {:.4f}'.format(epoch, loss/len(data)))
	return autoencoder


def plot_ae_components(data,autoencoder,device,label):
	for i, (x,y) in enumerate(data):
		z = autoencoder.encoder(x.to(device))
		z = z.to('cpu').detach().numpy()
		plt.scatter(z[:, 0], z[:, 1],c=y, cmap='tab10')
	plt.savefig("../output/"+label+".png");plt.close()

def get_encoded_z(df,autoencoder,device):
	x = df.iloc[:,1:-1].values.astype('float32')
	x = torch.from_numpy(x).to(device)
	z = autoencoder.encoder(x.to(device))
	df_z = pd.DataFrame(z.to('cpu').detach().numpy())
	df_z.columns = ["z"+str(i)for i in df_z.columns]
	return df_z

