from xml.dom.expatbuilder import theDOMImplementation
import torch; torch.manual_seed(0)
import torch.nn as nn
import pytorch_lightning as pl
from pytorch_lightning.plugins import DDPPlugin
from distribution import _dirichlet_multinomial as dirmult
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
		ss1 = self.fc1(xx1)

		xx2 = torch.log1p(xx2)
		ss2 = self.fc2(xx2)

		mm1 = self.z1_mean(ss1)
		lv1 = torch.clamp(self.z1_lnvar(ss1),-4.0,4.0)
		z1 = dirmult.reparameterize(mm1,lv1)

		mm2 = self.z2_mean(ss2)
		lv2 = torch.clamp(self.z2_lnvar(ss2),-4.0,4.0)
		z2 = dirmult.reparameterize(mm2,lv2)

		z = (z1 + z2) / 2

		return z,mm1,lv1,mm2,lv2

class ETMDecoder(nn.Module):
	def __init__(self,latent_dims,out_dims1, out_dims2,jitter=.1):
		super(ETMDecoder, self).__init__()

		self.p_beta_l_mean = nn.Parameter(torch.randn(latent_dims,out_dims1)*jitter)
		self.p_beta_r_mean = nn.Parameter(torch.randn(latent_dims,out_dims2)*jitter)

		self.p_beta_l_lnvar = nn.Parameter(torch.zeros(latent_dims,out_dims1))
		self.p_beta_r_lnvar = nn.Parameter(torch.zeros(latent_dims,out_dims2))

		self.p_beta_l_bias = nn.Parameter(torch.randn(1,out_dims1)*jitter)
		self.p_beta_r_bias = nn.Parameter(torch.randn(1,out_dims2)*jitter)

		self.lmbda_l = nn.Parameter(torch.randn(1)*jitter)
		self.lmbda_r = nn.Parameter(torch.randn(1)*jitter)

		self.l_smax = nn.LogSoftmax(dim=-1)

	def forward(self,zz):

		theta = torch.exp(self.l_smax(zz))

		z_beta_l, z_beta_r = self.get_beta()
		beta_l = z_beta_l.add(self.p_beta_l_bias)
		beta_r = z_beta_r.add(self.p_beta_r_bias)


		alpha_l = torch.exp(self.lmbda_l) * torch.exp(torch.clamp(torch.mm(theta,beta_l),-10,10))
		alpha_r = torch.exp(self.lmbda_r) * torch.exp(torch.clamp(torch.mm(theta,beta_r),-10,10))

		return self.p_beta_l_mean, self.p_beta_l_lnvar, self.p_beta_r_mean, self.p_beta_r_lnvar, theta, alpha_l, alpha_r
	
	def get_beta(self):

		lv_l = torch.clamp(self.p_beta_l_lnvar,-5.0,5.0) 
		lv_r = torch.clamp(self.p_beta_r_lnvar,-5.0,5.0)

		z_beta_l = dirmult.reparameterize(self.p_beta_l_mean,lv_l) 
		z_beta_r = dirmult.reparameterize(self.p_beta_r_mean,lv_r) 

		return z_beta_l, z_beta_r


class ETM(nn.Module):
	def __init__(self,input_dims1,input_dims2,latent_dims,layers1,layers2) :
		super(ETM,self).__init__()
		self.encoder = ETMEncoder(input_dims1,input_dims2,latent_dims,layers1,layers2)
		self.decoder = ETMDecoder(latent_dims, input_dims1, input_dims2)

	def forward(self,xx1,xx2):
		zz,m1,v1,m2,v2 = self.encoder(xx1,xx2)
		blm,blv,brm,brv,theta,alpha_l,alpha_r = self.decoder(zz)
		return m1,v1,m2,v2,blm,blv,brm,brv,theta,alpha_l,alpha_r

class LitETM(pl.LightningModule):

	def __init__(self,batch_size,input_dims1,input_dims2,latent_dims,layers1, layers2,lossf):
		super(LitETM,self).__init__()
		self.etm = ETM(input_dims1,input_dims2,latent_dims,layers1,layers2)
		self.data_size = 155e3 * 100
		self.lossf = lossf

	def forward(self,xx1,xx2):
		m1,v1,m2,v2,blm,blv,brm,brv,theta,alpha_l,alpha_r = self.etm(xx1,xx2)
		return m1,v1,m2,v2,blm,blv,brm,brv,theta,alpha_l,alpha_r

	def configure_optimizers(self):
		optimizer = torch.optim.Adam(self.parameters(), lr=0.01)
		return optimizer

	def training_step(self,batch):

		x1 , x2 = batch

		x1 = x1.reshape(x1.shape[0]*x1.shape[1],x1.shape[2])
		x2 = x2.reshape(x2.shape[0]*x2.shape[1],x2.shape[2])

		m1,v1,m2,v2,blm,blv,brm,brv,theta,alpha_l,alpha_r = self.etm(x1,x2)

		loglikloss_l = dirmult.log_likelihood(x1,alpha_l)
		loglikloss_r =  dirmult.log_likelihood(x2,alpha_r)
		loglikloss = loglikloss_l + loglikloss_r
		kl1 = dirmult.kl_loss(m1,v1)
		kl2 = dirmult.kl_loss(m2,v2)
		klb1 = dirmult.kl_loss(blm,blv)
		klb2 = dirmult.kl_loss(brm,brv)

		loss = torch.mean((kl1 + kl2)-loglikloss) + torch.sum(klb1)/self.data_size + torch.sum(klb2)/self.data_size

		f = open(self.lossf, 'a')
		f.write(str(torch.mean(loglikloss_l).to('cpu').item()) + '\t' +
				str(torch.mean(loglikloss_r).to('cpu').item()) + '\t' +
				str(torch.mean(loglikloss).to('cpu').item()) + '\t' +
		        str( torch.mean(kl1).to('cpu').item()) + '\t' + 
		        str( torch.mean(kl2).to('cpu').item()) + '\t' + 
				str((torch.sum(klb1)/self.data_size).to('cpu').item()) + '\t'+ 
				str((torch.sum(klb2)/self.data_size).to('cpu').item()) + '\n')
		f.close()

		return loss




