import torch; torch.manual_seed(0)
import torch.nn as nn
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
		indexes = self.index.get_nns_by_vector(vector.tolist(),k)
		return [self.labels[i][0] for i in indexes]

def get_NNmodels(df):
	model_list = {}
	for i in range(len(df['topic'].unique())):
		model_ann = ApproxNN(df.loc[df['topic']=='h'+str(i),['h'+str(x) for x in range(len(df['topic'].unique()))]].to_numpy(),df.loc[df['topic']=='h'+str(i),['cell']].values )
		model_ann.build()
		model_list[i] = model_ann
	return model_list

def get_neighbours(df,model_list,nbrsize):
	
	topic_len = len(df['topic'].unique())
	
	nbr_dict={}
	for idx in range(df.shape[0]):
		neighbours = np.array([ model_list[i].query(df.iloc[idx,2:].values,k=nbrsize) for i in range(topic_len)]).flatten()
		nbr_dict[idx] = df[df['cell'].isin(neighbours)].index.values

		if idx %1000 ==0:
			logger.info(idx)

	return nbr_dict

def generate_neighbours(args):

	args_home = os.environ['args_home']

	l_fname = args_home+args.input+args.raw_l_data
	r_fname = args_home+args.input+args.raw_r_data
	
	df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_h.tsv.gz',sep='\t',compression='gzip')

	df_l = pd.read_pickle(l_fname)
	df_r = pd.read_pickle(r_fname)
	df_l = df_l[df_l['index'].isin(df_h['cell'].values)]
	df_r = df_r[df_r['index'].isin(df_h['cell'].values)]

	dflatent = pd.merge(df_l['index'],df_h,how='left',left_on='index',right_on='cell')
	dflatent['cell'] = dflatent.iloc[:,2:].idxmax(axis=1)
	dflatent = dflatent.rename(columns={'cell':'topic','index':'cell'})

	model_list = get_NNmodels(dflatent)
	nbr_size = args.lr_model['train']['nbr_size']

	nbr_dict = get_neighbours(dflatent,model_list,nbr_size)

	pd.DataFrame(nbr_dict).T.to_pickle(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_nbr.pkl')

class SparseData():
	def __init__(self,indptr,indices,vals,shape,label,nbrmat,device):
		self.indptr = indptr
		self.indices = indices
		self.vals = vals
		self.shape = shape
		self.label = label
		self.nbrmat = nbrmat
		self.device = device

class SparseDataset(Dataset):
	def __init__(self, sparse_data):
		self.indptr = sparse_data.indptr
		self.indices = sparse_data.indices
		self.sparsemat = sparse_data.vals
		self.shape = sparse_data.shape
		self.device = sparse_data.device
		self.label = sparse_data.label
		self.nbrmat = sparse_data.nbrmat

	def __len__(self):
		return self.shape[0]

	def __getitem__(self, idx):

		self.device = torch.cuda.current_device()

		c_i = torch.zeros((self.shape[1],), dtype=torch.float32,requires_grad=False,device=self.device)
		ind1,ind2 = self.indptr[idx],self.indptr[idx+1]
		c_i[self.indices[ind1:ind2].long()] = self.sparsemat[ind1:ind2]
		ci_mat = c_i.expand(self.nbrmat.shape[1],1,c_i.shape[0]).squeeze(1)

		nbr_idxs = self.nbrmat[idx]
		cj_mat = []
		for i1,i2 in zip(self.indptr[nbr_idxs],self.indptr[nbr_idxs+1]):
			c_j = torch.zeros((self.shape[1],), dtype=torch.float32, requires_grad=False,device=self.device)
			c_j[self.indices[i1:i2].long()] = self.sparsemat[i1:i2]
			cj_mat.append(c_j)
		cj_mat = torch.stack(cj_mat).to(self.device)

		return ci_mat,cj_mat 


class load_data_lit(pl.LightningDataModule):

	def __init__(self,sparse_data,sparse_label,neighbour_data, batch_size,device):
		super().__init__()
		self.sparse_data = sparse_data
		self.sparse_label = sparse_label
		self.neighbour_data = neighbour_data
		self.batch_size = batch_size
		self.device = device

	def train_dataloader(self):

		npzarrs = np.load(self.sparse_data,allow_pickle=True)
		s = sparse.csr_matrix( (npzarrs['val'].astype(np.int32), (npzarrs['idx'], npzarrs['idy']) ),shape=npzarrs['shape'] )

		indptr = torch.tensor(s.indptr.astype(np.int32), dtype=torch.int32, device=self.device)
		indices = torch.tensor(s.indices.astype(np.int32), dtype=torch.int32, device=self.device)
		vals = torch.tensor(s.data.astype(np.int32), dtype=torch.int32, device=self.device)
		shape = tuple(npzarrs['shape'])

		metadat = np.load(self.sparse_label,allow_pickle=True)
		label = metadat['idx']

		df_nbr = pd.read_pickle(self.neighbour_data)

		nbrmat = torch.tensor(df_nbr.values.astype(np.compat.long),requires_grad=False).to(self.device)

		spdata = SparseData(indptr,indices,vals,shape,label,nbrmat,self.device)

		return DataLoader(SparseDataset(spdata), batch_size=self.batch_size, shuffle=True)

def load_data(sparse_data,sparse_label,neighbour_data, batch_size,device):

		npzarrs = np.load(sparse_data,allow_pickle=True)
		s = sparse.csr_matrix( (npzarrs['val'].astype(np.int32), (npzarrs['idx'], npzarrs['idy']) ),shape=npzarrs['shape'] )

		indptr = torch.tensor(s.indptr.astype(np.int32), dtype=torch.int32, device=device)
		indices = torch.tensor(s.indices.astype(np.int32), dtype=torch.int32, device=device)
		vals = torch.tensor(s.data.astype(np.int32), dtype=torch.float32, device=device)
		shape = tuple(npzarrs['shape'])

		metadat = np.load(sparse_label,allow_pickle=True)
		label = metadat['idx']

		df_nbr = pd.read_pickle(neighbour_data)

		nbrmat = torch.tensor(df_nbr.values.astype(np.compat.long),requires_grad=False).to(device)

		spdata = SparseData(indptr,indices,vals,shape,label,nbrmat,device)

		return DataLoader(SparseDataset(spdata), batch_size=batch_size, shuffle=True)

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

		self.p_alpha= nn.Parameter(torch.randn(latent_dims,out_dims1)*jitter)
		self.alpha_bias= nn.Parameter(torch.randn(1,out_dims1)*jitter)

		self.p_beta= nn.Parameter(torch.randn(latent_dims,out_dims1,out_dims2)*jitter)
		self.beta_bias= nn.Parameter(torch.randn(1,out_dims2)*jitter)

	def forward(self,zz,xx1):

		hh = self.l_smax(zz)

		alpha = self.l_smax(self.alpha_bias.add(self.p_alpha))
		px = torch.mm(torch.exp(hh),torch.exp(alpha))
			
		beta = self.l_smax(self.beta_bias.add(self.p_beta))
		beta_x = torch.matmul(xx1,torch.exp(beta)).sum(1)
		py = torch.mm(torch.exp(hh),beta_x)

		return px,py,hh

class ETM(nn.Module):
	def __init__(self,input_dims1,input_dims2,latent_dims,layers1, layers2) :
		super(ETM,self).__init__()
		self.encoder = ETMEncoder(input_dims1,input_dims2,latent_dims,layers1,layers2)
		self.decoder = ETMDecoder(latent_dims, input_dims1, input_dims2)

	def forward(self,xx1,xx2):
		zz,m1,v1,m2,v2 = self.encoder(xx1,xx2)
		px,py,h = self.decoder(zz,xx1)
		return px,py,m1,v1,m2,v2,h

class LitETM(pl.LightningModule):

	def __init__(self,input_dims1,input_dims2,latent_dims,layers1, layers2,lossf):
		super(LitETM,self).__init__()
		self.etm = ETM(input_dims1,input_dims2,latent_dims,layers1,layers2)
		self.lossf = lossf

	def forward(self,xx1,xx2):
		px,py,m1,v1,m2,v2,h = self.etm(xx1,xx2)
		return px,py,m1,v1,m2,v2,h

	def configure_optimizers(self):
		optimizer = torch.optim.Adam(self.parameters(), lr=0.01)
		return optimizer

	def training_step(self,batch):

		x1 , x2 = batch

		x1 = x1.squeeze(0)
		x2 = x2.squeeze(0)

		px,py,m1,v1,m2,v2,h = self.etm(x1,x2)

		x = torch.cat((x1,x2),1)
		loglikloss_x = etm_llik(x1,px)
		loglikloss_y = etm_llik(x2,py)
		loglikloss = loglikloss_x + loglikloss_y
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

	nbr_size = args.lr_model['train']['nbr_size']
	batch_size = args.lr_model['train']['batch_size']
	l_rate = args.lr_model['train']['l_rate']
	epochs = args.lr_model['train']['epochs']
	layers1 = args.lr_model['train']['layers1']
	layers2 = args.lr_model['train']['layers2']
	latent_dims = args.lr_model['train']['latent_dims']
	input_dims1 = args.lr_model['train']['input_dims1']
	input_dims2 = args.lr_model['train']['input_dims2']


	args_home = os.environ['args_home']
	device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

	sparse_data = args_home+args.input+args.nbr_model['sparse_data']
	sparse_label = args_home+args.input+args.nbr_model['sparse_label']
	neighbour_data = args_home+ args.output+ args.lr_model['out']+args.nbr_model['mfile']+'_nbr.pkl'

	dl = load_data(sparse_data,sparse_label,neighbour_data, batch_size, device)

	train_dataloader =  dl.train_dataloader()

	logging.info('Input dimension - ligand is '+ str(input_dims1))
	logging.info('Input dimension - receptor is '+ str(input_dims2))
	model = LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,model_file+'loss.txt')
	logging.info(model)

	trainer = pl.Trainer(
	max_epochs=epochs,
	accelerator='gpu',
	plugins= DDPPlugin(find_unused_parameters=False),
	gradient_clip_val=0.5,
	progress_bar_refresh_rate=50,
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

	model = LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,'temp.txt')

	model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']
	model.load_state_dict(torch.load(model_file+'_ietm.torch'))
	model.eval()

	beta1,beta1_bias,beta2,beta2_bias =  None,None,None,None
	
	for n,p in model.named_parameters():
		print(n)
		if n == 'etm.decoder.lbeta1':
			beta1=p
		elif n == 'etm.decoder.lbeta2':
			beta2=p
		elif n == 'etm.decoder.lbeta1_bias':
			beta1_bias=p
		elif n == 'etm.decoder.lbeta2_bias':
			beta2_bias=p
		

	beta_smax = nn.LogSoftmax(dim=-1)
	beta1 = torch.exp(beta_smax(beta1))
	beta2 = torch.exp(beta_smax(beta2))

	df_beta1 = pd.DataFrame(beta1.to('cpu').detach().numpy())
	df_beta1.to_csv(model_file+'_ietm_beta1.tsv.gz',sep='\t',index=False,compression='gzip')

	df_beta2 = pd.DataFrame(beta2.to('cpu').detach().numpy())
	df_beta2.to_csv(model_file+'_ietm_beta2.tsv.gz',sep='\t',index=False,compression='gzip')
	
	df_beta1_bias = pd.DataFrame(beta1_bias.to('cpu').detach().numpy())
	df_beta1_bias.to_csv(model_file+'_ietm_beta1_bias.tsv.gz',sep='\t',index=False,compression='gzip')
	
	df_beta2_bias = pd.DataFrame(beta2_bias.to('cpu').detach().numpy())
	df_beta2_bias.to_csv(model_file+'_ietm_beta2_bias.tsv.gz',sep='\t',index=False,compression='gzip')

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
		r = 0
		for x1,x2 in data:
			x1 = x1.squeeze(0)
			x2 = x2.squeeze(0)
			opt.zero_grad()
			px,py,m1,v1,m2,v2,h = etm(x1,x2)

			loglikloss_x = etm_llik(x1,px)
			loglikloss_y = etm_llik(x2,py)
			loglikloss = loglikloss_x + loglikloss_y
			kl1 = kl_loss(m1,v1)
			kl2 = kl_loss(m2,v2)
			kl = kl1+kl2

			train_loss = torch.mean(kl-loglikloss).to("cpu")
			train_loss.backward()
			opt.step()
			loss += train_loss.item()
			r += 1
			print(epoch,r)
		logger.info('====> Epoch: {} Average loss: {:.4f}'.format(epoch, loss/len(data)))
		
		loss_values.append(loss/len(data))
	
	return loss_values


def run_model(args,model_file):

	nbr_size = args.lr_model['train']['nbr_size']
	batch_size = args.lr_model['train']['batch_size']
	l_rate = args.lr_model['train']['l_rate']
	epochs = args.lr_model['train']['epochs']
	layers1 = args.lr_model['train']['layers1']
	layers2 = args.lr_model['train']['layers2']
	latent_dims = args.lr_model['train']['latent_dims']
	input_dims1 = args.lr_model['train']['input_dims1']
	input_dims2 = args.lr_model['train']['input_dims2']


	args_home = os.environ['args_home']
	device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

	sparse_data = args_home+args.input+args.nbr_model['sparse_data']
	sparse_label = args_home+args.input+args.nbr_model['sparse_label']
	neighbour_data = args_home+ args.output+ args.lr_model['out']+args.nbr_model['mfile']+'_nbr.pkl'

	dl = load_data(sparse_data,sparse_label,neighbour_data, batch_size, device)

	model = ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
	logging.info(model)
	loss_values = train(model,dl,device,epochs,l_rate)

	torch.save(model.state_dict(), model_file+"etm.torch")
	dflv = pd.DataFrame(loss_values)
	dflv.to_csv(model_file+"loss.txt",index=False)




