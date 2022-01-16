

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
		output = F.relu(self.linear1(x))
		z = F.relu(self.linear2(output))
		return z

class VariationalEncoder(nn.Module):
	def __init__(self, input_dims, latent_dims):
		super(VariationalEncoder, self).__init__()
		self.linear1 = nn.Linear(input_dims, 128)
		self.linear2 = nn.Linear(128, latent_dims)
		self.linear3 = nn.Linear(128, latent_dims)

		## TODO need to look up pytorch documentation for cuda sampling call
		self.N = torch.distributions.Normal(0, 1)
		self.N.loc = self.N.loc.cuda() 
		self.N.scale = self.N.scale.cuda()
		self.kl = 0

	def forward(self, x):
		output = F.relu(self.linear1(x))
		mu =  self.linear2(output)
		## TODO do we need to pass sigma from relu as well if exp is used
		sigma = torch.exp(self.linear3(output))
		z = mu + sigma*self.N.sample(mu.shape)
		self.kl = (sigma**2 + mu**2 - torch.log(sigma) - 1/2).sum()
		return z


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
		self.encoder = Encoder(input_dims,latent_dims)
		self.decoder = Decoder(input_dims,latent_dims)
	def forward(self, x):
		z = self.encoder(x)
		return self.decoder(z)

class VariationalAutoencoder(nn.Module):
	def __init__(self,input_dims,latent_dims):
		super(VariationalAutoencoder, self).__init__()
		self.encoder = VariationalEncoder(input_dims,latent_dims)
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


def vaetrain(autoencoder, data,device, epochs=50):
	opt = torch.optim.Adam(autoencoder.parameters())
	mse = nn.MSELoss()
	for epoch in range(epochs):
		loss = 0
		for indx,(x,y) in enumerate(data):
			x = x.to(device)
			opt.zero_grad()
			x_hat = autoencoder(x)
			train_loss = mse(x_hat,x) + autoencoder.encoder.kl
			train_loss.backward()
			opt.step()
			loss += train_loss.item()
		if epoch % 5 == 0:  
			print('====> Epoch: {} Average loss: {:.4f}'.format(epoch, loss/len(data)))
	return autoencoder

