import torch; torch.manual_seed(0)
import torch.nn as nn
from distribution import _multinomial as st
import logging
logger = logging.getLogger(__name__)

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

class ETMEncoder(nn.Module):
	def __init__(self,input_dims,latent_dims,layers):
		super(ETMEncoder, self).__init__()
		self.fc = Stacklayers(input_dims,layers)
		self.z_mean = nn.Linear(latent_dims,latent_dims)
		self.z_lnvar = nn.Linear(latent_dims,latent_dims)

	def forward(self, xx):

		xx = torch.log1p(xx)
		xx = xx/torch.sum(xx,dim=-1,keepdim=True)
		ss = self.fc(xx)

		mm = self.z_mean(ss)
		lv = torch.clamp(self.z_lnvar(ss),-4.0,4.0)
		z = st.reparameterize(mm,lv)
		return z,mm,lv

class ETMDecoder(nn.Module):
	def __init__(self,latent_dims,out_dims,jitter=.1):
		super(ETMDecoder, self).__init__()
		self.lbeta_bias= nn.Parameter(torch.randn(1,out_dims)*jitter)
		self.lbeta = nn.Parameter(torch.randn(latent_dims,out_dims)*jitter)
		self.beta = nn.LogSoftmax(dim=-1)
		self.hid = nn.LogSoftmax(dim=-1)

	def forward(self, zz):
		beta = self.beta(self.lbeta_bias.add(self.lbeta))
		hh = self.hid(zz)
		return torch.mm(torch.exp(hh),torch.exp(beta)), hh

class ETM(nn.Module):
	def __init__(self,input_dims,latent_dims,layers):
		super(ETM,self).__init__()
		self.encoder = ETMEncoder(input_dims,latent_dims,layers)
		self.decoder = ETMDecoder(latent_dims, input_dims)

	def forward(self,xx):
		zz,m,v = self.encoder(xx)
		pr,h = self.decoder(zz)
		return pr,m,v,h

def train(etm,data,epochs,l_rate):
	logger.info('Starting training....')
	opt = torch.optim.Adam(etm.parameters(),lr=l_rate)
	loss_values = []
	loss_values_sep = []
	for epoch in range(epochs):
		loss = 0
		loss_ll = 0
		loss_kl = 0
		for x,y in data:
			opt.zero_grad()
			recon,m,v,h = etm(x)
			loglikloss = st.log_likelihood(x,recon)
			kl = st.kl_loss(m,v)

			ll_l = torch.mean(loglikloss).to('cpu')
			kl_ml = torch.mean(kl).to('cpu')

			train_loss = torch.mean(kl-loglikloss).to('cpu')
			train_loss.backward()

			opt.step()

			loss += train_loss.item()
			loss_ll += ll_l.item()
			loss_kl += kl_ml.item()

		if epoch % 10 == 0:
			logger.info('====> Epoch: {} Average loss: {:.4f}'.format(epoch, loss/len(data)))

		loss_values.append(loss/len(data))
		loss_values_sep.append((loss_ll/len(data),loss_kl/len(data)))


	return loss_values,loss_values_sep

