import torch; torch.manual_seed(0)
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
from scipy import sparse
import numpy as np
import pandas as pd
import os
import logging
logger = logging.getLogger(__name__)


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

class TabularDataset(Dataset):
	def __init__(self, x,y):
		self.X = x
		self.y = y
	def __len__(self):
		return len(self.X)
	def __getitem__(self, idx):
		return [self.X[idx], self.y[idx]]

def load_data(x1,x2,device,bath_size):
	x1 = x1.astype('float32')
	x1 = torch.from_numpy(x1).to(device)
	x2 = x2.astype('float32')
	x2 = torch.from_numpy(x2).to(device)
	return DataLoader(TabularDataset(x1,x2), batch_size=bath_size, shuffle=True)

class SparseData():
	def __init__(self,indptr,indices,data,shape,immuneindx,otherindx,label):
		self.indptr = indptr
		self.indices = indices
		self.data = data
		self.shape = shape
		self.immuneindx = immuneindx
		self.otherindx = otherindx
		self.label = label

class SparseTabularDataset(Dataset):
	def __init__(self, sparse_data,device):
		self.indptr = sparse_data.indptr
		self.indices = sparse_data.indices
		self.data = sparse_data.data
		self.shape = sparse_data.shape
		self.immuneindx = sparse_data.immuneindx
		self.otherindx = sparse_data.otherindx
		self.label = sparse_data.label
		self.device = device

	def __len__(self):
		return self.shape[0]

	def __getitem__(self, idx):

		cell = torch.zeros((self.shape[1],), dtype=torch.int32, device=self.device)
		ind1,ind2 = self.indptr[idx],self.indptr[idx+1]
		cell[self.indices[ind1:ind2].long()] = self.data[ind1:ind2]

		cell_i = cell[self.immuneindx.long()]
		cell_o = cell[self.otherindx.long()]

		return cell_i,cell_o,self.label[idx]

def load_sparse_data(data_file,meta_file,immuneindex_file,device,bath_size):

	logger.info('loading sparse data...\n'+
		data_file +'\n'+
		meta_file +'\n' +
		immuneindex_file )

	npzarrs = np.load(data_file,allow_pickle=True)
	s = sparse.csr_matrix( (npzarrs['val'].astype(np.int32), (npzarrs['idx'], npzarrs['idy']) ),shape=npzarrs['shape'] )

	device = torch.device(device)
	indptr = torch.tensor(s.indptr.astype(np.int32), dtype=torch.int32, device=device)
	indices = torch.tensor(s.indices.astype(np.int32), dtype=torch.int32, device=device)
	val = torch.tensor(s.data.astype(np.int32), dtype=torch.int32, device=device)
	shape = tuple(npzarrs['shape'])

	metadat = np.load(meta_file,allow_pickle=True)
	label = metadat['idx']

	immunedat = np.load(immuneindex_file,allow_pickle=True)
	immuneindx = torch.tensor(np.array([x[0] for x in immunedat['idx']]).astype(np.int32), dtype=torch.int32, device=device)

	otherindx = torch.tensor(np.array([x for x in range(shape[1]) if x not in immuneindx]).astype(np.int32),dtype=torch.int32, device=device)

	spdata = SparseData(indptr,indices,val,shape,immuneindx,otherindx,label)

	return DataLoader(SparseTabularDataset(spdata,device), batch_size=bath_size, shuffle=True)

def train(etm,data,epochs,l_rate):
	logger.info('Starting training....')
	opt = torch.optim.Adam(etm.parameters(),lr=l_rate)
	loss_values = []
	loss_values_sep = []
	for epoch in range(epochs):
		loss = 0
		loss_ll = 0
		loss_kl1 = 0
		loss_kl2 = 0
		for x1,x2,y in data:
			opt.zero_grad()
			recon,m1,v1,m2,v2,h = etm(x1,x2)
			x = torch.cat((x1,x2),1)
			loglikloss = etm_llik(x,recon)
			kl1 = kl_loss(m1,v1)
			kl2 = kl_loss(m2,v2)

			ll_l = torch.mean(loglikloss).to('cpu')
			kl_ml1 = torch.mean(kl1).to('cpu')
			kl_ml2 = torch.mean(kl2).to('cpu')

			train_loss = torch.mean((kl1 + kl2)-loglikloss).to('cpu')
			train_loss.backward()

			opt.step()

			loss += train_loss.item()
			loss_ll += ll_l.item()
			loss_kl1 += kl_ml1.item()
			loss_kl2 += kl_ml2.item()

		if epoch % 10 == 0:
			logger.info('====> Epoch: {} Average loss: {:.4f}'.format(epoch, loss/len(data)))

		loss_values.append(loss/len(data))
		loss_values_sep.append((loss_ll/len(data),loss_kl1/len(data),
		loss_kl2/len(data)))


	return loss_values,loss_values_sep

def get_latent(data,model,model_file):

	for xx1,xx2,label in data: break
	zz,m1,v1,m2,v2 = model.encoder(xx1,xx2)
	pr,hh = model.decoder(zz)
	hh = torch.exp(hh)

	df_z = pd.DataFrame(zz.to('cpu').detach().numpy())
	df_z.columns = ['zz'+str(i)for i in df_z.columns]
	df_z['cell'] = label
	df_z = df_z[ ['cell']+[x for x in df_z.columns if x not in['cell','sample']]]
	df_z.to_csv(model_file+'etm_zz_data.tsv',sep='\t',index=False,compression='gzip')

	df_h = pd.DataFrame(hh.to('cpu').detach().numpy())
	df_h.columns = ['hh'+str(i)for i in df_h.columns]
	df_h['cell'] = label
	df_h = df_h[ ['cell']+[x for x in df_h.columns if x not in['cell','sample']]]
	df_h.to_csv(model_file+'etm_hh_data.tsv',sep='\t',index=False,compression='gzip')

	beta1 =  None
	beta2 =  None
	for n,p in model.named_parameters():
		if n == 'decoder.lbeta1':
			beta1=p
		elif n == 'decoder.lbeta2':
			beta2=p
	beta_smax = nn.LogSoftmax(dim=-1)
	beta1 = torch.exp(beta_smax(beta1))
	beta2 = torch.exp(beta_smax(beta2))

	df_beta1 = pd.DataFrame(beta1.to('cpu').detach().numpy())
	df_beta1.to_csv(model_file+'etm_beta1_data.tsv',sep='\t',index=False,compression='gzip')

	df_beta2 = pd.DataFrame(beta2.to('cpu').detach().numpy())
	df_beta2.to_csv(model_file+'etm_beta2_data.tsv',sep='\t',index=False,compression='gzip')

def run_model(args):

	args_home = os.environ['args_home']
	logger.info('Starting model training...')

	data_file = args_home+args.input+args.nbr_model['sparse_data']
	meta_file = args_home+args.input+args.nbr_model['sparse_label']
	immuneindex_file = args_home+args.input+args.nbr_model['sparse_immune_label']

	batch_size = args.nbr_model['train']['batch_size']
	l_rate = args.nbr_model['train']['l_rate']
	epochs = args.nbr_model['train']['epochs']
	layers1 = args.nbr_model['train']['layers1']
	layers2 = args.nbr_model['train']['layers2']
	latent_dims = args.nbr_model['train']['latent_dims']

	device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
	data = load_sparse_data(data_file,meta_file,immuneindex_file, device,batch_size)

	input_dims1 = len(data.dataset.immuneindx)
	input_dims2 = data.dataset.shape[1] - input_dims1
	logging.info('Input dimension - immune is '+ str(input_dims1))
	logging.info('Input dimension - others is '+ str(input_dims2))

	model = ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
	logging.info(model)
	loss_values = train(model,data,epochs,l_rate)

	torch.save(model.state_dict(), model_file+'etm.torch')
	dflv = pd.DataFrame(loss_values[0])
	dflv.to_csv(model_file+'loss.txt',index=False)
	dflv = pd.DataFrame(loss_values[1])
	dflv.to_csv(model_file+'loss2.txt',index=False)

def eval_model(args):

	logging.info('Starting model inference...')

	data_file = args_home+args.input+args.nbr_model['sparse_data']
	meta_file = args_home+args.input+args.nbr_model['sparse_label']
	immuneindex_file = args_home+args.input+args.nbr_model['sparse_immune_label']

	batch_size = args.nbr_model['eval']['batch_size']
	layers1 = args.nbr_model['train']['layers1']
	layers2 = args.nbr_model['train']['layers2']
	latent_dims = args.nbr_model['train']['latent_dims']

	model_file = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']

	device = 'cpu'

	data = load_sparse_data(data_file,meta_file,immuneindex_file, device,batch_size)

	print('data size is -',data.dataset.shape)

	input_dims1 = len(data.dataset.immuneindx)
	input_dims2 = data.dataset.shape[1] - input_dims1

	model = ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
	model.load_state_dict(torch.load(model_file+'etm.torch'))
	model.eval()

	dfloss = pd.read_csv(model_file+'loss.txt')

	get_latent(data,model,model_file)

