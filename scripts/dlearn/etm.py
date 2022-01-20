
import torch; torch.manual_seed(0)
import torch.nn as nn
import torch.nn.functional as F


class Stacklayers(nn.Module):
	def __init__(self,input_size,layers):
		super(Stacklayers, self).__init__()
		self.layers = nn.ModuleList()
		self.input_size = input_size
		for next_size,activfunc in layers:
			self.layers.append(nn.Linear(input_size,next_size))
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
		elif act == 'softplus':
			act = nn.Softplus()
		elif act == 'rrelu':
			act = nn.RReLU()
		elif act == 'leakyrelu':
			act = nn.LeakyReLU()
		elif act == 'elu':
			act = nn.ELU()
		elif act == 'selu':
			act = nn.SELU()
		elif act == 'glu':
			act = nn.GLU()
		else:
			print('Defaulting to tanh activations...')
			act = nn.Tanh()
		return act 


class ETMEncoder(nn.Module):
	def __init__(self,input_dims,latent_dims,layers):
		super(ETMEncoder, self).__init__()
		self.fc = Stacklayers(input_dims,layers)
		hidden_dims = layers[len(layers)][0]
		self.z_mu = nn.Linear(hidden_dims,latent_dims)
		self.z_var = nn.Linear(hidden_dims,latent_dims)
		self.latent_dim = latent_dims

	def forward(self, xx):
		# xx = ## get normalized of xx
		ss = self.fc(torch.log1p(xx))
		mm = self.z_mu(ss)
		lv = torch.exp(self.z_var(ss))
		z = self.normal_sample(mm,lv)
		return [z,mm,lv]
	
	def normal_sample(self,mean,logvar):
		std = torch.exp(logvar/2.)
		eps = torch.randn_like(std)
		return eps.mul_(std).add_(mean)

	def kl_loss(self,mean,logvar):
		return  -0.5 * torch.sum(1. + logvar - mean.pow(2) - logvar.exp(), dim=-1).mean()
	

class ETMDecoder(nn.Module):
	def __init__(self,out_dims,latent_dims,jitter=.1):
		super(ETMDecoder, self).__init__()
		self.lbeta = nn.parameter(torch.randn(latent_dims,out_dims)*jitter)
		self.beta = nn.LogSoftmax(2)
		self.hid = nn.LogSoftmax(2)

	def forward(self, zz):
		beta = self.beta(self.lbeta)
		hh = self.hid(zz)
		return torch.mm(torch.exp(hh),torch.exp(beta))

class ETM(nn.Module):
	def __init__(self,input_dims,out_dims,latent_dims,layers) :
		super(ETM,self).__init__()
		self.encoder = ETMEncoder(input_dims,latent_dims,layers)
		self.decoder = ETMDecoder(out_dims,latent_dims)
		
	def forward(self,xx):
		zz = self.encoder(xx)
		recon = self.decoder(zz)
		return [recon,zz[1],zz[2]]
		
mlnn = Stacklayers(1500, \
[(128,"relu"),\
(100,"relu"),\
(50,"relu"),\
(32,"relu")])
print(mlnn)