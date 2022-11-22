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
		
		self.beta_bias= nn.Parameter(torch.randn(1,out_dims)*jitter)
		self.beta_mean = nn.Parameter(torch.randn(latent_dims,out_dims)*jitter)
		self.beta_lnvar = nn.Parameter(torch.zeros(latent_dims,out_dims))

		self.lsmax = nn.LogSoftmax(dim=-1)

	def forward(self, zz):
		
		theta = torch.exp(self.lsmax(zz))
		
		z_beta = self.get_beta()
		beta = z_beta.add(self.beta_bias)

		return self.beta_mean,self.beta_lnvar,theta,beta

	def get_beta(self):
		lv = torch.clamp(self.beta_lnvar,-5.0,5.0)
		z_beta = st.reparameterize(self.beta_mean,lv) 
		return z_beta

class ETM(nn.Module):
	def __init__(self,input_dims,latent_dims,layers):
		super(ETM,self).__init__()
		self.encoder = ETMEncoder(input_dims,latent_dims,layers)
		self.decoder = ETMDecoder(latent_dims, input_dims)

	def forward(self,xx):
		zz,m,v = self.encoder(xx)
		bm,bv,theta,beta = self.decoder(zz)
		return m,v,bm,bv,theta,beta

def train(etm,data,epochs,l_rate,batch_size):
	logger.info('Starting training....')
	opt = torch.optim.Adam(etm.parameters(),lr=l_rate)
	loss_values = []
	loss_values_sep = []
	data_size = data.dataset.shape[0]
	for epoch in range(epochs):
		loss = 0
		loss_ll = 0
		loss_kl = 0
		loss_klb = 0
		for x,y in data:
			opt.zero_grad()
			m,v,bm,bv,theta,beta = etm(x)
			alpha = torch.exp(torch.clamp(torch.mm(theta,beta),-10,10))
			loglikloss = st.multi_dir_log_likelihood(x,alpha)
			kl = st.kl_loss(m,v)
			klb = st.kl_loss(bm,bv)

			train_loss = torch.mean(kl -loglikloss).add(torch.sum(klb)/data_size)
			train_loss.backward()

			opt.step()

			ll_l = torch.mean(loglikloss).to('cpu')
			kl_l = torch.mean(kl).to('cpu')
			klb_l = torch.sum(klb).to('cpu')

			loss += train_loss.item()
			loss_ll += ll_l.item()
			loss_kl += kl_l.item()
			loss_klb += klb_l.item()/data_size

		if epoch % 10 == 0:
			logger.info('====> Epoch: {} Average loss: {:.4f}'.format(epoch, loss/len(data)))

		loss_values.append(loss/len(data))
		loss_values_sep.append((loss_ll/len(data),loss_kl/len(data),loss_klb/len(data)))


	return loss_values,loss_values_sep

