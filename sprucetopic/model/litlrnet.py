import torch; torch.manual_seed(0)
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
import pytorch_lightning as pl
from pytorch_lightning.plugins import DDPPlugin
import numpy as np
import pandas as pd
import random
import annoy
import os
import logging
logger = logging.getLogger(__name__)

class ApproxNN():
	def __init__(self, data, labels):
		self.dimension = data.shape[1]
		self.data = data.astype('float32')
		self.labels = labels

	def build(self, number_of_trees=50):
		self.index = annoy.AnnoyIndex(self.dimension,'angular')
		for i, vec in enumerate(self.data):
			self.index.add_item(i, vec.tolist())
		self.index.build(number_of_trees)

	def query(self, vector, k):
		return self.index.get_nns_by_vector(vector.tolist(),k)

class LRDataset(Dataset):
	def __init__(self,lmat,rmat,Alr,nbrsize,hmat,model_ann) :
		self.lmat = lmat
		self.rmat = rmat
		self.lrmat = Alr
		self.model_ann = model_ann
		self.nbrsize = nbrsize
		self.dflatent = hmat

	def __len__(self):
		return self.dflatent.shape[0]

	def __getitem__(self, idx):

		cm_l = self.lmat[idx].unsqueeze(0)
		cm_r = self.rmat[idx].unsqueeze(0)

		cm_lr = torch.mm(cm_l,self.lrmat).mul(cm_r)
		cm_rl = torch.mm(cm_r,torch.t(self.lrmat)).mul(cm_l)

		neighbours = self.model_ann.query(self.dflatent.iloc[idx,2:].values,k=45)
		# nbr_idxs = self.dflatent.iloc[neighbours].groupby('topic').head(1).index

		cn_l = self.lmat[neighbours]
		cn_r = self.rmat[neighbours]

		return 	cm_lr + torch.mm(cn_l,self.lrmat).mul(cn_r) , cm_rl + torch.mm(cn_r,torch.t(self.lrmat)).mul(cn_l)

class load_data(pl.LightningDataModule):

	def __init__(self, args,nbr_size, batch_size):
		super().__init__()
		self.args = args
		self.nbr_size = nbr_size
		self.batch_size = batch_size

	def train_dataloader(self):

		device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

		args_home = os.environ['args_home']

		l_fname = args_home+self.args.input+self.args.raw_l_data
		r_fname = args_home+self.args.input+self.args.raw_r_data
		lr_fname = args_home+self.args.input+self.args.raw_lr_data


		df_h = pd.read_csv(args_home+self.args.output+self.args.nbr_model["out"]+self.args.nbr_model["mfile"]+"etm_hh_data.tsv",sep="\t",compression='gzip')

		df_l = pd.read_pickle(l_fname)
		df_r = pd.read_pickle(r_fname)
		df_l = df_l[df_l["index"].isin(df_h["cell"].values)]
		df_r = df_r[df_r["index"].isin(df_h["cell"].values)]

		dflatent = pd.merge(df_l['index'],df_h,how='left',left_on='index',right_on='cell')
		dflatent['cell'] = dflatent.iloc[:,2:].idxmax(axis=1)
		dflatent = dflatent.rename(columns={'cell':'topic'})

		model_ann = ApproxNN(dflatent.iloc[:,2:].to_numpy(), dflatent.iloc[:,0].values)
		model_ann.build()

		# h_mat = torch.tensor(dfjoin.iloc[:,2:].values.astype(np.float32),requires_grad=False).to(device)
		lmat = torch.tensor(df_l.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(device)
		rmat = torch.tensor(df_r.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(device)

		df_lr = pd.read_pickle(lr_fname)
		df_lr = df_lr.loc[df_l.columns[1:],df_r.columns[1:]]
		Alr = torch.tensor(df_lr.values.astype(np.float32),requires_grad=False).to(device)

		return DataLoader(LRDataset(lmat,rmat,Alr,self.nbr_size,dflatent,model_ann), batch_size=self.batch_size, shuffle=True)

def reparameterize(mean,lnvar):
	sig = torch.exp(lnvar/2.)
	eps = torch.randn_like(sig)
	return eps.mul_(sig).add_(mean)

def etm_llik(xx,pr, eps=1e-8):
	return torch.sum(xx * torch.log(pr+eps),dim=-1)

def kl_loss(mean,lnvar):
	return  -0.5 * torch.sum(1. + lnvar - torch.pow(mean,2) - torch.exp(lnvar), dim=-1)

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
		xx1 = xx1/torch.sum(xx1,dim=-1,keepdim=True)
		xx2 = xx2/torch.sum(xx2,dim=-1,keepdim=True)

		ss1 = self.fc1(torch.log1p(xx1))
		ss2 = self.fc2(torch.log1p(xx2))

		mm1 = self.z1_mean(ss1)
		lv1 = torch.clamp(self.z1_lnvar(ss1),-4.0,4.0)
		z1 = reparameterize(mm1,lv1)

		mm2 = self.z2_mean(ss2)
		lv2 = torch.clamp(self.z2_lnvar(ss2),-4.0,4.0)
		z2 = reparameterize(mm2,lv2)

		z = (z1 + z2) / 2

		return z,mm1,lv1,mm2,lv2


class ETMDecoder(nn.Module):
	def __init__(self,latent_dims,out_dims1, out_dims2,jitter=.1):
		super(ETMDecoder, self).__init__()
		self.lbeta1= nn.Parameter(torch.randn(latent_dims,out_dims1)*jitter)
		self.lbeta2= nn.Parameter(torch.randn(latent_dims,out_dims2)*jitter)
		self.beta = nn.LogSoftmax(dim=-1)
		self.hid = nn.LogSoftmax(dim=-1)

	def forward(self, zz):
		beta1 = self.beta(self.lbeta1)
		beta2 = self.beta(self.lbeta2)

		hh = self.hid(zz)

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

	def __init__(self,input_dims1,input_dims2,latent_dims,layers1, layers2):
		super(LitETM,self).__init__()
		self.etm = ETM(input_dims1,input_dims2,latent_dims,layers1,layers2)

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
		loglikloss = etm_llik(x,recon)
		kl1 = kl_loss(m1,v1)
		kl2 = kl_loss(m2,v2)
		loss = torch.mean((kl1 + kl2)-loglikloss)

		self.log("train_loss", loss)

		return loss

def run_model(args,):

	args_home = os.environ['args_home']


	nbr_size = args.lr_model['train']['nbr_size']
	batch_size = args.lr_model['train']['batch_size']
	l_rate = args.lr_model['train']['l_rate']
	epochs = args.lr_model['train']['epochs']
	layers1 = args.lr_model['train']['layers1']
	layers2 = args.lr_model['train']['layers2']
	latent_dims = args.lr_model['train']['latent_dims']


	dl = load_data(args,nbr_size,batch_size)

	train_dataloader =  dl.train_dataloader()

	input_dims1 = 539
	input_dims2 = 498
	logging.info('Input dimension - ligand is '+ str(input_dims1))
	logging.info('Input dimension - receptor is '+ str(input_dims2))
	model = LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2)
	logging.info(model)

	trainer = pl.Trainer(
	max_epochs=epochs,
	accelerator='gpu',
	plugins= DDPPlugin(find_unused_parameters=False),
	precision=16,
	progress_bar_refresh_rate=10)
	
	trainer.fit(model,train_dataloader)

# if __name__ == '__main__':
# 	mode=sys.argv[1]
# 	run_model(mode)


