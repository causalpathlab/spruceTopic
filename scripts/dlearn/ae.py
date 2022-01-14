

## in progress..
## compiled from pytorch tutorials..
## https://avandekleut.github.io/vae/
## https://github.com/lschmiddey/Autoencoder


import torch; torch.manual_seed(0)
import torchvision
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.data import Dataset


class Encoder(nn.Module):
	def __init__(self, input_dims, latent_dims):
		super(Encoder, self).__init__()
		self.linear1 = nn.Linear(input_dims, 128)
		self.linear2 = nn.Linear(128, latent_dims)
	def forward(self, x):
		x = torch.flatten(x, start_dim=1)
		x = F.relu(self.linear1(x))
		return self.linear2(x)

class Decoder(nn.Module):
	def __init__(self, input_dims, latent_dims):
		super(Decoder, self).__init__()
		self.linear1 = nn.Linear(latent_dims, 128)
		self.linear2 = nn.Linear(128, input_dims)
	def forward(self, z):
		z = F.relu(self.linear1(z))
		z = torch.sigmoid(self.linear2(z))
		return z

class Autoencoder(nn.Module):
	def __init__(self,input_dims,latent_dims):
		super(Autoencoder, self).__init__()
		self.encoder = Encoder(input_dims,latent_dims)
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
	# opt = torch.optim.Adam(autoencoder.parameters())
	opt = torch.optim.SGD(autoencoder.parameters(),lr=0.01)
	# loss = nn.MSELoss(reduction="sum")
	train_loss = []
	for epoch in range(epochs):
		for indx,(x,y) in enumerate(data):
			x = x.to(device)
			opt.zero_grad()
			x_hat = autoencoder(x)
			loss = ((x - x_hat)**2).sum()
			train_loss.append(loss.item())
			loss.backward()
			opt.step()
		if epoch % 30 == 0:        
			print('====> Epoch: {} Average loss: {:.4f}'.format(epoch, train_loss[-1]))
	return autoencoder


