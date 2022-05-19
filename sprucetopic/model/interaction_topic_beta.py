import torch.nn as nn
import torch; torch.manual_seed(0)
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
import pytorch_lightning as pl
from pytorch_lightning.plugins import DDPPlugin
import numpy as np
import pandas as pd
from scipy import sparse
import annoy
import os
import logging
logger = logging.getLogger(__name__)

class load_data_lit(pl.LightningDataModule):

	def __init__(self,f_latent_h,f_l,f_r,f_neighbour, batch_size,device):
		super().__init__()
		self.f_latent_h = f_latent_h
		self.f_l = f_l
		self.f_r = f_r
		self.f_neighbour = f_neighbour
		self.batch_size = batch_size
		self.device = device

	def train_dataloader(self):
		df_h = pd.read_csv(self.f_latent_h,sep='\t')
		df_l = pd.read_pickle(self.f_l)
		df_r = pd.read_pickle(self.f_r)
		df_l = df_l[df_l['index'].isin(df_h['cell'].values)]
		df_r = df_r[df_r['index'].isin(df_h['cell'].values)]

		df_nbr = pd.read_pickle(self.f_neighbour)
		nbrmat = torch.tensor(df_nbr.values.astype(np.compat.long),requires_grad=False).to(self.device)

		lmat = torch.tensor(df_l.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(self.device)
		rmat = torch.tensor(df_r.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(self.device)

		return DataLoader(LRDataset(lmat,rmat,nbrmat), batch_size=self.batch_size, shuffle=True)


def load_data(f_latent_h,f_l,f_r,f_neighbour, batch_size,device):
	df_h = pd.read_csv(f_latent_h,sep='\t')
	df_l = pd.read_pickle(f_l)
	df_r = pd.read_pickle(f_r)
	df_l = df_l[df_l['index'].isin(df_h['cell'].values)]
	df_r = df_r[df_r['index'].isin(df_h['cell'].values)]

	df_nbr = pd.read_pickle(f_neighbour)
	nbrmat = torch.tensor(df_nbr.values.astype(np.compat.long),requires_grad=False).to(device)

	lmat = torch.tensor(df_l.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(device)
	rmat = torch.tensor(df_r.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(device)

	return DataLoader(LRDataset(lmat,rmat,nbrmat), batch_size=batch_size, shuffle=True)

class LRDataset(Dataset):
	def __init__(self,lmat,rmat,nbrmat) :
		self.lmat = lmat
		self.rmat = rmat
		self.nbrmat = nbrmat

	def __len__(self):
		return len(self.nbrmat)

	def __getitem__(self, idx):

		c_i = self.lmat[idx].unsqueeze(0)
		ci_mat = c_i.expand(self.nbrmat.shape[1],1,c_i.shape[1]).squeeze(1)
		
		cj_mat = self.rmat[self.nbrmat[idx]]

		return ci_mat,cj_mat 

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

		xx1 = torch.log1p(xx1)
		# xx1 = xx1/torch.sum(xx1,dim=-1,keepdim=True)
		ss1 = self.fc1(xx1)

		xx2 = torch.log1p(xx2)
		# xx2 = xx2/torch.sum(xx2,dim=-1,keepdim=True)
		ss2 = self.fc2(xx2)

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

		self.l_smax = nn.LogSoftmax(dim=-1)
		self.p_beta= nn.Parameter(torch.randn(latent_dims,out_dims1,out_dims2)*jitter)
		self.beta_bias= nn.Parameter(torch.randn(1,out_dims2)*jitter)

	def forward(self,zz,xx1):
		hh = self.l_smax(zz)			
		beta = self.l_smax(self.beta_bias.add(self.p_beta))
		beta_x = torch.matmul(xx1,torch.exp(beta)).sum(1)
		py = torch.mm(torch.exp(hh),beta_x)

		return py,hh

class ETM(nn.Module):
	def __init__(self,input_dims1,input_dims2,latent_dims,layers1, layers2) :
		super(ETM,self).__init__()
		self.encoder = ETMEncoder(input_dims1,input_dims2,latent_dims,layers1,layers2)
		self.decoder = ETMDecoder(latent_dims, input_dims1, input_dims2)

	def forward(self,xx1,xx2):
		zz,m1,v1,m2,v2 = self.encoder(xx1,xx2)
		py,h = self.decoder(zz,xx1)
		return py,m1,v1,m2,v2,h

class LitETM(pl.LightningModule):

	def __init__(self,input_dims1,input_dims2,latent_dims,layers1, layers2,lossf):
		super(LitETM,self).__init__()
		self.etm = ETM(input_dims1,input_dims2,latent_dims,layers1,layers2)
		self.lossf = lossf

	def forward(self,xx1,xx2):
		py,m1,v1,m2,v2,h = self.etm(xx1,xx2)
		return py,m1,v1,m2,v2,h

	def configure_optimizers(self):
		optimizer = torch.optim.Adam(self.parameters(), lr=0.01)
		return optimizer

	def training_step(self,batch):

		x1 , x2 = batch

		x1 = x1.reshape(x1.shape[0]*x1.shape[1],x1.shape[2])
		x2 = x2.reshape(x2.shape[0]*x2.shape[1],x2.shape[2])

		py,m1,v1,m2,v2,h = self.etm(x1,x2)

		loglikloss = etm_llik(x2,py)
		kl1 = kl_loss(m1,v1)
		kl2 = kl_loss(m2,v2)
		loss = torch.mean((kl1 + kl2)-loglikloss)

		self.log('train_loss', loss)

		ll_l = torch.mean(loglikloss).to('cpu').item()
		kl_ml1 = torch.mean(kl1).to('cpu').item()
		kl_ml2 = torch.mean(kl2).to('cpu').item()

		f = open(self.lossf, 'a')
		f.write(str(ll_l) + ';' + str(kl_ml1) + ';'+ str(kl_ml2) + '\n')
		f.close()

		return loss

def run_model_lit(args,model_file):

	batch_size = args.lr_model['train']['batch_size']
	l_rate = args.lr_model['train']['l_rate']
	epochs = args.lr_model['train']['epochs']
	layers1 = args.lr_model['train']['layers1']
	layers2 = args.lr_model['train']['layers2']
	latent_dims = args.lr_model['train']['latent_dims']
	input_dims1 = args.lr_model['train']['input_dims1']
	input_dims2 = args.lr_model['train']['input_dims2']


	args_home = os.environ['args_home']

	f_l= args_home+args.input+args.raw_l_data
	f_r = args_home+args.input+args.raw_r_data
	f_latent_h = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_h.tsv.gz'
	f_neighbour = args_home+args.output+args.lr_model['out']+args.nbr_model['mfile']+'_nbr.pkl'

	device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

	dl = load_data_lit(f_latent_h,f_l,f_r,f_neighbour, batch_size,device)

	train_dataloader =  dl.train_dataloader()

	logging.info('Input dimension - ligand is '+ str(input_dims1))
	logging.info('Input dimension - receptor is '+ str(input_dims2))
	model = LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,model_file+'loss.txt')
	logging.info(model)
	for n,p in model.named_parameters():
		logging.info(n)
	trainer = pl.Trainer(
	max_epochs=epochs,
	accelerator='gpu',
	plugins= DDPPlugin(find_unused_parameters=False),
	gradient_clip_val=0.5,
	progress_bar_refresh_rate=500,
	enable_checkpointing=False)
	
	trainer.fit(model,train_dataloader)

	torch.save(model.state_dict(), model_file+'_ietm.torch')

def load_model(args):

	args_home = os.environ['args_home']

	layers1 = args.lr_model['train']['layers1']
	layers2 = args.lr_model['train']['layers2']
	latent_dims = args.lr_model['train']['latent_dims']
	input_dims1 = args.lr_model['train']['input_dims1']
	input_dims2 = args.lr_model['train']['input_dims2']

	# model = LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,'temp.txt')
	model = ETM(input_dims1,input_dims2,latent_dims,layers1,layers2)

	model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']
	model.load_state_dict(torch.load(model_file+'_ietm.torch'))
	model.eval()

	alpha,alpha_bias,beta,beta_bias =  None,None,None,None
	
	for n,p in model.named_parameters():
		print(n)
		if n == 'decoder.p_alpha':
			alpha=p
		elif n == 'decoder.p_beta':
			beta=p
		elif n == 'decoder.alpha_bias':
			alpha_bias=p
		elif n == 'decoder.beta_bias':
			beta_bias=p
		

	beta_smax = nn.LogSoftmax(dim=-1)
	alpha = torch.exp(beta_smax(alpha))
	beta = torch.exp(beta_smax(beta))

	df_alpha = pd.DataFrame(alpha.to('cpu').detach().numpy())
	df_alpha.to_csv(model_file+'_ietm_alpha.tsv.gz',sep='\t',index=False,compression='gzip')

	df_beta = pd.DataFrame(beta.to('cpu').detach().numpy())
	df_beta.to_csv(model_file+'_ietm_beta.tsv.gz',sep='\t',index=False,compression='gzip')
	
	df_alpha_bias = pd.DataFrame(alpha_bias.to('cpu').detach().numpy())
	df_alpha_bias.to_csv(model_file+'_ietm_alpha_bias.tsv.gz',sep='\t',index=False,compression='gzip')
	
	df_beta_bias = pd.DataFrame(beta_bias.to('cpu').detach().numpy())
	df_beta_bias.to_csv(model_file+'_ietm_beta_bias.tsv.gz',sep='\t',index=False,compression='gzip')

def eval_model(args):

	args_home = os.environ['args_home']

	layers1 = args.lr_model['train']['layers1']
	layers2 = args.lr_model['train']['layers2']
	latent_dims = args.lr_model['train']['latent_dims']
	input_dims1 = args.lr_model['train']['input_dims1']
	input_dims2 = args.lr_model['train']['input_dims2']

	model = LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,'_txt_')

	model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']
	model.load_state_dict(torch.load(model_file+'_ietm.torch'))
	model.eval()


	l_fname = args_home+args.input+args.raw_l_data
	r_fname = args_home+args.input+args.raw_r_data
	lr_fname = args_home+args.input+args.raw_lr_data

	df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_h.tsv.gz',sep='\t',compression='gzip')

	df_l = pd.read_pickle(l_fname)
	df_r = pd.read_pickle(r_fname)
	df_l = df_l[df_l['index'].isin(df_h['cell'].values)]
	df_r = df_r[df_r['index'].isin(df_h['cell'].values)]
	df_lr = pd.read_pickle(lr_fname)
	df_lr = df_lr.loc[df_l.columns[1:],df_r.columns[1:]]
	df_nbr = pd.read_pickle(args_home+args.output+args.lr_model['out']+args.nbr_model['mfile']+'_nbr.pkl')

	# device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
	device='cpu'
	nbrmat = torch.tensor(df_nbr.values.astype(np.compat.long),requires_grad=False).to(device)
	lmat = torch.tensor(df_l.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(device)
	rmat = torch.tensor(df_r.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(device)
	lrmat = torch.tensor(df_lr.values.astype(np.float32),requires_grad=False).to(device)

	topics = []
	for idx in range(df_nbr.shape[0]):

		cm_l = lmat[idx].unsqueeze(0)
		cm_r = rmat[idx].unsqueeze(0)

		cm_lr = torch.mm(cm_l,lrmat).mul(cm_r)
		cm_rl = torch.mm(cm_r,torch.t(lrmat)).mul(cm_l)

		nbr_idxs = nbrmat[idx]
		
		cn_l = lmat[nbr_idxs]
		cn_r = rmat[nbr_idxs]

		lprime,rprime =  cm_lr + torch.mm(cn_l,lrmat).mul(cn_r) , cm_rl + torch.mm(cn_r,torch.t(lrmat)).mul(cn_l)

		pr,m1,v1,m2,v2,h = model.etm(lprime,rprime)

		h_smax = nn.LogSoftmax(dim=-1)
		h = torch.exp(h_smax(h))

		topics.append(list(pd.DataFrame(h.detach().numpy()).idxmax(axis=1).values))
	
	df_it = pd.DataFrame(topics)
	df_it['cell'] = df_l['index']
	df_it = df_it[['cell']+[x for x in df_it.columns[:-1]]]
	df_it.to_csv(model_file+'_ietm_interaction_states.tsv.gz',sep='\t',index=False,compression='gzip')



def train(etm, data, device, epochs,l_rate):
	logger.info("Starting training....")
	opt = torch.optim.Adam(etm.parameters(),lr=l_rate)
	loss_values = []
	for epoch in range(epochs):
		loss = 0
		for x1,x2 in data:

			x1 = x1.reshape(x1.shape[0]*x1.shape[1],x1.shape[2])
			x2 = x2.reshape(x2.shape[0]*x2.shape[1],x2.shape[2])

			opt.zero_grad()
			py,m1,v1,m2,v2,h = etm(x1,x2)

			loglikloss_y = etm_llik(x2,py)
			loglikloss = loglikloss_y
			kl1 = kl_loss(m1,v1)
			kl2 = kl_loss(m2,v2)
			kl = kl1+kl2

			train_loss = torch.mean(kl-loglikloss).to("cpu")
			train_loss.backward()
			opt.step()
			loss += train_loss.item()
		logger.info('====> Epoch: {} Average loss: {:.4f}'.format(epoch, loss/len(data)))
		
		loss_values.append(loss/len(data))
	
	return loss_values


def run_model(args,model_file):

	batch_size = args.lr_model['train']['batch_size']
	l_rate = args.lr_model['train']['l_rate']
	epochs = args.lr_model['train']['epochs']
	layers1 = args.lr_model['train']['layers1']
	layers2 = args.lr_model['train']['layers2']
	latent_dims = args.lr_model['train']['latent_dims']
	input_dims1 = args.lr_model['train']['input_dims1']
	input_dims2 = args.lr_model['train']['input_dims2']


	args_home = os.environ['args_home']

	f_l= args_home+args.input+args.raw_l_data
	f_r = args_home+args.input+args.raw_r_data
	f_latent_h = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_h.tsv.gz'
	f_neighbour = args_home+args.output+args.lr_model['out']+args.nbr_model['mfile']+'_nbr.pkl'

	device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
	dl = load_data(f_latent_h,f_l,f_r,f_neighbour, batch_size,device)

	model = ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
	logging.info(model)
	loss_values = train(model,dl,device,epochs,l_rate)

	torch.save(model.state_dict(), model_file+"etm.torch")
	dflv = pd.DataFrame(loss_values)
	dflv.to_csv(model_file+"loss.txt",index=False)