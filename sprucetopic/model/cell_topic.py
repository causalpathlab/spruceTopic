import torch; torch.manual_seed(0)
import torch.nn as nn
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
		z = reparameterize(mm,lv)
		return z,mm,lv

# class ETMDecoder(nn.Module):
# 	def __init__(self,latent_dims,out_dims,jitter=.1):
# 		super(ETMDecoder, self).__init__()
# 		self.lbeta = nn.Parameter(torch.randn(latent_dims,out_dims)*jitter)
# 		self.beta = nn.LogSoftmax(dim=-1)
# 		self.hid = nn.LogSoftmax(dim=-1)

# 	def forward(self, zz):
# 		beta = self.beta(self.lbeta)
# 		hh = self.hid(zz)
# 		return torch.mm(torch.exp(hh),torch.exp(beta)), hh
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

class SparseData():
	def __init__(self,indptr,indices,vals,shape,label):
		self.indptr = indptr
		self.indices = indices
		self.vals = vals
		self.shape = shape
		self.label = label

class SparseTabularDataset(Dataset):
	def __init__(self, sparse_data,device):
		self.indptr = sparse_data.indptr
		self.indices = sparse_data.indices
		self.vals = sparse_data.vals
		self.shape = sparse_data.shape
		self.device = device
		self.label = sparse_data.label

	def __len__(self):
		return self.shape[0]

	def __getitem__(self, idx):

		cell = torch.zeros((self.shape[1],), dtype=torch.int32, device=self.device)
		ind1,ind2 = self.indptr[idx],self.indptr[idx+1]
		cell[self.indices[ind1:ind2].long()] = self.vals[ind1:ind2]

		return cell,self.label[idx]

def load_sparse_data(data_file,meta_file,device,bath_size):

	logger.info('loading sparse data...\n'+
		data_file +'\n'+ meta_file )

	npzarrs = np.load(data_file,allow_pickle=True)
	s = sparse.csr_matrix( (npzarrs['val'].astype(np.int32), (npzarrs['idx'], npzarrs['idy']) ),shape=npzarrs['shape'] )

	device = torch.device(device)
	indptr = torch.tensor(s.indptr.astype(np.int32), dtype=torch.int32, device=device)
	indices = torch.tensor(s.indices.astype(np.int32), dtype=torch.int32, device=device)
	vals = torch.tensor(s.data.astype(np.int32), dtype=torch.int32, device=device)
	shape = tuple(npzarrs['shape'])

	metadat = np.load(meta_file,allow_pickle=True)
	label = metadat['idx']

	spdata = SparseData(indptr,indices,vals,shape,label)

	return DataLoader(SparseTabularDataset(spdata,device), batch_size=bath_size, shuffle=True)

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
			loglikloss = etm_llik(x,recon)
			kl = kl_loss(m,v)

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


def run_model(args,model_file):

	args_home = os.environ['args_home']
	logger.info('Starting model training...')

	data_file = args_home+args.input+args.nbr_model['sparse_data']
	meta_file = args_home+args.input+args.nbr_model['sparse_label']

	batch_size = args.nbr_model['train']['batch_size']
	l_rate = args.nbr_model['train']['l_rate']
	epochs = args.nbr_model['train']['epochs']
	layers = args.nbr_model['train']['layers']
	latent_dims = args.nbr_model['train']['latent_dims']

	device = args.nbr_model['train']['device']
	
	# torch.device('cuda' if torch.cuda.is_available() else 'cpu')
	data = load_sparse_data(data_file,meta_file, device,batch_size)

	input_dims = data.dataset.shape[1]
	logging.info('Input dimension is '+ str(input_dims))

	model = ETM(input_dims,latent_dims,layers).to(device)
	logging.info(model)
	loss_values = train(model,data,epochs,l_rate)

	torch.save(model.state_dict(), model_file+'_netm.torch')
	dflv = pd.DataFrame(loss_values[0])
	dflv.to_csv(model_file+'_netm_loss.txt.gz',index=False,compression='gzip')
	dflv = pd.DataFrame(loss_values[1])
	dflv.to_csv(model_file+'_netm_loss2.txt.gz',index=False,compression='gzip')

def eval_model(args):

	logging.info('Starting model inference...')
	args_home = os.environ['args_home']

	data_file = args_home+args.input+args.nbr_model['sparse_data']
	label_file = args_home+args.input+args.nbr_model['sparse_label']

	batch_size = args.nbr_model['eval']['batch_size']
	layers = args.nbr_model['train']['layers']
	latent_dims = args.nbr_model['train']['latent_dims']

	model_file = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']

	device = 'cpu'

	data = load_sparse_data(data_file,label_file,device,batch_size)

	print('data size is -',data.dataset.shape)

	input_dims = data.dataset.shape[1]

	model = ETM(input_dims,latent_dims,layers).to(device)
	model.load_state_dict(torch.load(model_file+'_netm.torch'))
	model.eval()

	get_latent(data,model,model_file,label_file)

def get_latent(data,model,model_file,label_file):

	for xx,y in data: break
	zz,m,v = model.encoder(xx)
	pr,hh = model.decoder(zz)
	hh = torch.exp(hh)

	df_z = pd.DataFrame(zz.to('cpu').detach().numpy())
	df_z.columns = ['z'+str(i)for i in df_z.columns]
	df_z['cell'] = y
	df_z = df_z[ ['cell']+[x for x in df_z.columns if x not in['cell']]]
	df_z.to_csv(model_file+'_netm_z.tsv.gz',sep='\t',index=False,compression='gzip')

	df_h = pd.DataFrame(hh.to('cpu').detach().numpy())
	df_h.columns = ['h'+str(i)for i in df_h.columns]
	df_h['cell'] = y
	df_h = df_h[ ['cell']+[x for x in df_h.columns if x not in['cell']]]
	df_h.to_csv(model_file+'_netm_h.tsv.gz',sep='\t',index=False,compression='gzip')

	beta =  None
	for n,p in model.named_parameters():
		if n == 'decoder.lbeta':
			beta=p
	beta_smax = nn.LogSoftmax(dim=-1)
	beta = torch.exp(beta_smax(beta))

	df_beta1 = pd.DataFrame(beta.to('cpu').detach().numpy())
	df_beta1.to_csv(model_file+'_netm_beta.tsv.gz',sep='\t',index=False,compression='gzip')
