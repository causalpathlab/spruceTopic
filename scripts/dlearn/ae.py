
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
	def __init__(self,input_size,layers):
		super(Stacklayers, self).__init__()
		self.layers = nn.ModuleList()
		self.input_size = input_size
		for next_l in layers:
			self.layers.append(nn.Linear(self.input_size,next_l))
			self.layers.append(self.get_activation())
			nn.BatchNorm1d(next_l)
			self.input_size = next_l
		
	def forward(self, input_data):
		for layer in self.layers:
			input_data = layer(input_data)
		return input_data

	def get_activation(self):
		return nn.ReLU()


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
	def __init__(self,input_dims,out_dims,latent_dims,layers):
		super(Autoencoder, self).__init__()
		self.encoder = Stacklayers(input_dims,layers)
		self.decoder = Decoder(out_dims,latent_dims)
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

def load_data(x,y,device,batch_size):
	x = x.astype('float32')
	x = torch.from_numpy(x).to(device)
	return DataLoader(TabularDataset(x,y), batch_size=batch_size, shuffle=True)

def train(autoencoder, data,device, epochs,lr):
	opt = torch.optim.Adam(autoencoder.parameters(),lr=lr)
	mse = nn.MSELoss()
	loss_values =[]
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
		loss_values.append(loss/len(data))
	return loss_values


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


def get_encoded_h(df,model,device,title,loss_values):
	import pandas as pd
	import matplotlib.pylab as plt
	import seaborn as sns
	from matplotlib.cm import ScalarMappable

	x = df.iloc[:,1:-1].values.astype('float32')
	x = torch.from_numpy(x).to(device)
	
	z = model.encoder(x.to(device))

	df_z = pd.DataFrame(z.to('cpu').detach().numpy())
	df_z.columns = ["z"+str(i)for i in df_z.columns]
	
	data_color = range(len(df_z.columns))
	data_color = [x / max(data_color) for x in data_color] 

	custom_map = plt.cm.get_cmap('coolwarm') #one of the color schemas stored
	custom = custom_map(data_color)  #mapping the color info to the variable custom
	df_z.plot(kind='bar', stacked=True, color=custom,figsize=(25,10))
	plt.ylabel("hidden state proportion", fontsize=18)
	plt.xlabel("samples", fontsize=22)
	plt.xticks([])
	plt.title(title,fontsize=25)
	plt.savefig("../output/ae_z_"+title+".png");plt.close()

	plt.plot(loss_values)
	plt.ylabel("loss", fontsize=18)
	plt.xlabel("epochs", fontsize=22)
	plt.title(title,fontsize=25)
	plt.savefig("../output/ae_z_"+title+"_loss.png");plt.close()
	# df_z.to_csv("../output/hh_"+title+"_data.csv",sep="\t",index=False)


