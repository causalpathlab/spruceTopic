import torch; torch.manual_seed(0)
import torch.nn as nn
import pytorch_lightning as pl
from pytorch_lightning.plugins import DDPPlugin
from distribution import _multinomial as st
import numpy as np
import pandas as pd
import os
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
	def __init__(self,input_dims1, input_dims2,latent_dims,layers1,layers2):
		super(ETMEncoder, self).__init__()
		self.fc1 = Stacklayers(input_dims1,layers1)
		self.fc2 = Stacklayers(input_dims2,layers2)
		self.z1_mean = nn.Linear(latent_dims,latent_dims)
		self.z1_lnvar = nn.Linear(latent_dims,latent_dims)
		self.z2_mean = nn.Linear(latent_dims,latent_dims)
		self.z2_lnvar = nn.Linear(latent_dims,latent_dims)

	def forward(self, xx1, xx2):

		xx1 = torch.log1p(xx1)
		# xx1 = xx1/torch.sum(xx1,dim=-1,keepdim=True)
		ss1 = self.fc1(xx1)

		xx2 = torch.log1p(xx2)
		# xx2 = xx2/torch.sum(xx2,dim=-1,keepdim=True)
		ss2 = self.fc2(xx2)

		mm1 = self.z1_mean(ss1)
		lv1 = torch.clamp(self.z1_lnvar(ss1),-4.0,4.0)
		z1 = st.reparameterize(mm1,lv1)

		mm2 = self.z2_mean(ss2)
		lv2 = torch.clamp(self.z2_lnvar(ss2),-4.0,4.0)
		z2 = st.reparameterize(mm2,lv2)

		z = (z1 + z2) / 2

		return z,mm1,lv1,mm2,lv2

class ETMDecoder(nn.Module):
	def __init__(self,latent_dims,out_dims1, out_dims2,jitter=.1):
		super(ETMDecoder, self).__init__()

		self.p_beta1= nn.Parameter(torch.randn(latent_dims,out_dims1)*jitter)
		self.p_beta2= nn.Parameter(torch.randn(latent_dims,out_dims2)*jitter)

		self.p_beta1_bias= nn.Parameter(torch.randn(1,out_dims1)*jitter)
		self.p_beta2_bias= nn.Parameter(torch.randn(1,out_dims2)*jitter)

		self.l_smax = nn.LogSoftmax(dim=-1)

	def forward(self, zz):
		beta1 = self.l_smax(self.p_beta1_bias.add(self.p_beta1))
		beta2 = self.l_smax(self.p_beta2_bias.add(self.p_beta2))

		hh = self.l_smax(zz)

		beta = torch.cat((beta1,beta2),1)

		return torch.mm(torch.exp(hh),torch.exp(beta)), hh

class ETM(nn.Module):
	def __init__(self,input_dims1,input_dims2,latent_dims,layers1, layers2) :
		super(ETM,self).__init__()
		self.encoder = ETMEncoder(input_dims1,input_dims2,latent_dims,layers1,layers2)
		self.decoder = ETMDecoder(latent_dims, input_dims1, input_dims2)

	def forward(self,xx1,xx2):
		zz,m1,v1,m2,v2 = self.encoder(xx1,xx2)
		pr,h = self.decoder(zz)
		return pr,m1,v1,m2,v2,h

class LitETM(pl.LightningModule):

	def __init__(self,input_dims1,input_dims2,latent_dims,layers1, layers2,lossf):
		super(LitETM,self).__init__()
		self.etm = ETM(input_dims1,input_dims2,latent_dims,layers1,layers2)
		self.lossf = lossf

	def forward(self,xx1,xx2):
		pr,m1,v1,m2,v2,h = self.etm(xx1,xx2)
		return pr,m1,v1,m2,v2,h

	def configure_optimizers(self):
		optimizer = torch.optim.Adam(self.parameters(), lr=0.01)
		return optimizer

	def training_step(self,batch):

		x1 , x2 = batch

		x1 = x1.reshape(x1.shape[0]*x1.shape[1],x1.shape[2])
		x2 = x2.reshape(x2.shape[0]*x2.shape[1],x2.shape[2])

		recon,m1,v1,m2,v2,h = self.etm(x1,x2)

		x = torch.cat((x1,x2),1)
		loglikloss = st.log_likelihood(x,recon)
		kl1 = st.kl_loss(m1,v1)
		kl2 = st.kl_loss(m2,v2)
		loss = torch.mean((kl1 + kl2)-loglikloss)

		self.log('train_loss', loss)

		ll_l = torch.mean(loglikloss).to('cpu').item()
		kl_ml1 = torch.mean(kl1).to('cpu').item()
		kl_ml2 = torch.mean(kl2).to('cpu').item()

		f = open(self.lossf, 'a')
		f.write(str(ll_l) + ';' + str(kl_ml1) + ';'+ str(kl_ml2) + '\n')
		f.close()

		return loss






