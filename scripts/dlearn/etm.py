
from cmath import log
import torch; torch.manual_seed(0)
import torch.nn as nn
import torch.nn.functional as F
import torchvision
from torch.utils.data import DataLoader
from torch.utils.data import Dataset


class Stacklayers(nn.Module):
	def __init__(self,input_size,layers):
		super(Stacklayers, self).__init__()
		self.layers = nn.ModuleList()
		self.input_size = input_size
		for next_size in layers:
			self.layers.append(nn.Linear(input_size,next_size))
			input_size = next_size
			self.layers.append(self.get_activation())
		
	def forward(self, input_data):
		for layer in self.layers:
			input_data = layer(input_data)
		return input_data

	def get_activation(self):
		return nn.ReLU()

class ETMEncoder(nn.Module):
	def __init__(self,input_dims,latent_dims,layers):
		super(ETMEncoder, self).__init__()
		self.fc = Stacklayers(input_dims,layers)
		hidden_dims = layers[len(layers)-1]
		self.latent_dim = latent_dims
		self.z_mu = nn.Linear(hidden_dims,latent_dims)
		self.z_var = nn.Linear(hidden_dims,latent_dims)

	def normal_sample(self,mean,logvar):
		std = torch.exp(0.5*logvar)
		eps = torch.randn_like(std)
		return eps.mul_(std).add_(mean)


	def forward(self, xx):
		## TODO normalization ? q_theta = qtheta/torch.sum(xx,dim-1)
		q_theta = self.fc(torch.log1p(xx))
		mu_theta = self.z_mu(q_theta)
		logsigma_theta = self.z_var(q_theta)

		zz = self.normal_sample(mu_theta,logsigma_theta)

		return zz, mu_theta,logsigma_theta


class ETMDecoder(nn.Module):
	def __init__(self,out_dims,latent_dims,jitter=.1):
		super(ETMDecoder, self).__init__()
		self.q_beta = nn.Parameter(torch.randn(latent_dims,out_dims)*jitter)

	def forward(self, zz):
		beta = F.softmax(self.q_beta,dim=-1)
		hh = F.softmax(zz,dim=-1)
		return torch.mm(torch.exp(hh),torch.exp(beta))

class ETM(nn.Module):
	def __init__(self,input_dims,out_dims,latent_dims,layers) :
		super(ETM,self).__init__()
		self.encoder = ETMEncoder(input_dims,latent_dims,layers)
		self.decoder = ETMDecoder(out_dims,latent_dims)
		
	def forward(self,xx):
		zz,mu_theta,logsigma_theta = self.encoder(xx)
		recon = self.decoder(zz)
		return recon,mu_theta,logsigma_theta

def etm_llik(xx,pr, eps=1e-8):
	# return torch.sum(xx* torch.log(pr+eps),dim=-1)
	recon_loss = -(pr * xx).sum(1)
	return recon_loss.mean()

def kl_loss(mean,logvar):
	return  -0.5 * torch.sum(1. + logvar - mean.pow(2) - logvar.exp(), dim=-1).mean()

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

def train(etm, data,device, epochs=300):
	opt = torch.optim.Adam(etm.parameters())
	for epoch in range(epochs):
		loss = 0
		for indx,(x,y) in enumerate(data):
			x = x.to(device)
			opt.zero_grad()
			x_hat,mu_theta,logsigma_theta = etm(x)
			loglikloss = etm_llik(x,x_hat)
			klloss = kl_loss(mu_theta,logsigma_theta)
			train_loss = loglikloss + klloss
			train_loss.backward()
			opt.step()
			loss += train_loss.item()
		if epoch % 30 == 0:  
			print('====> Epoch: {} Average loss: {:.4f}'.format(epoch, loss/len(data)))
	return etm
