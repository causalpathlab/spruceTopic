import torch; torch.manual_seed(0)
import torch.nn as nn
import torch.nn.functional as F
import torchvision
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
from scipy import sparse
import numpy as np 
import preprocess


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
		hidden_dims = layers[len(layers)-1]
		self.latent_dim = latent_dims
		self.z_mean = nn.Linear(hidden_dims,latent_dims)
		self.z_lnvar = nn.Linear(hidden_dims,latent_dims)

	def forward(self, xx):
		xx = xx/torch.sum(xx,dim=-1,keepdim=True)
		ss = self.fc(torch.log1p(xx))
		mm = self.z_mean(ss)
		lv = torch.clamp(self.z_lnvar(ss),-4.0,4.0)
		z = reparameterize(mm,lv)
		return z, mm,lv

class ETMDecoder(nn.Module):
	def __init__(self,out_dims,latent_dims,jitter=.1):
		super(ETMDecoder, self).__init__()
		self.lbeta = nn.Parameter(torch.randn(latent_dims,out_dims)*jitter)
		self.beta = nn.LogSoftmax(dim=-1)
		self.hid = nn.LogSoftmax(dim=-1)

	def forward(self, zz):
		beta = self.beta(self.lbeta)
		hh = self.hid(zz)
		return torch.mm(torch.exp(hh),torch.exp(beta)), torch.exp(hh)
	
class ETM(nn.Module):
	def __init__(self,input_dims,out_dims,latent_dims,layers) :
		super(ETM,self).__init__()
		self.encoder = ETMEncoder(input_dims,latent_dims,layers)
		self.decoder = ETMDecoder(out_dims,latent_dims)
		
	def forward(self,xx):
		zz,mean,lv = self.encoder(xx)
		pr,hh = self.decoder(zz)
		return pr,mean,lv,hh
	
class TabularDataset(Dataset):
	def __init__(self, x,y):
		self.X = x
		self.y = y	
	def __len__(self):
		return len(self.X)
	def __getitem__(self, idx):
		return [self.X[idx], self.y[idx]]

def load_data(x,y,device,bath_size):
	x = x.astype('float32')
	x = torch.from_numpy(x).to(device)
	return DataLoader(TabularDataset(x,y), batch_size=bath_size, shuffle=True)

class SparseData():
	def __init__(self,indptr,indices,data,shape,label):
		self.indptr = indptr
		self.indices = indices
		self.data = data
		self.shape = shape
		self.label = label

class SparseTabularDataset(Dataset):
	def __init__(self, sparse_data, device):
		self.indptr = sparse_data.indptr
		self.indices = sparse_data.indices
		self.data = sparse_data.data
		self.shape = sparse_data.shape
		self.label = sparse_data.label
		self.device = device

	def __len__(self):
		return self.shape[0]

	def __getitem__(self, idx):
		cell = torch.zeros((self.shape[1],), dtype=torch.int32, device=self.device)
		ind1,ind2 = self.indptr[idx],self.indptr[idx+1]
		cell[self.indices[ind1:ind2].long()] = self.data[ind1:ind2]
		return cell,self.label[idx]
	
	def dmat(self):
		s = torch.sparse_csr_tensor(self.indptr,self.indices,self.data,size=self.shape)
		return s.to_dense()

def load_sparse_data(data_file,meta_file,device,bath_size):

	npzarrs = np.load(data_file,allow_pickle=True)
	s = sparse.csr_matrix( (npzarrs['val'].astype(np.int32), (npzarrs['idx'], npzarrs['idy']) ),shape=npzarrs['shape'] )
	
	device = torch.device(device)
	indptr = torch.tensor(s.indptr.astype(np.int32), dtype=torch.int32, device=device)
	indices = torch.tensor(s.indices.astype(np.int32), dtype=torch.int32, device=device)
	val = torch.tensor(s.data.astype(np.int32), dtype=torch.int32, device=device)
	shape = tuple(npzarrs['shape'])

	metadat = np.load(meta_file,allow_pickle=True)
	label = metadat["idx"]
	spdata = SparseData(indptr,indices,val,shape,label)
	return DataLoader(SparseTabularDataset(spdata,device), batch_size=bath_size, shuffle=True)


def train(etm, data,device, epochs,l_rate):
	opt = torch.optim.Adam(etm.parameters(),lr=l_rate)
	loss_values = []
	for epoch in range(epochs):
		loss = 0
		for indx,(x,y) in enumerate(data):
			x = x.to(device)
			opt.zero_grad()
			recon,mean,lnvar,hh = etm(x)
			loglikloss = etm_llik(x,recon)
			kl = kl_loss(mean,lnvar)
			train_loss = torch.mean(kl-loglikloss).to("cpu")
			train_loss.backward()
			opt.step()
			loss += train_loss.item()
		if epoch % 100 == 0:  
			print('====> Epoch: {} Average loss: {:.4f}'.format(epoch, loss/len(data)))
		
		loss_values.append(loss/len(data))
	
	return loss_values


def get_encoded_h(x,label,model,device,title,loss_values):
	import pandas as pd
	import matplotlib.pylab as plt
	import seaborn as sns
	from matplotlib.cm import ScalarMappable
	
	zz,mean,lv = model.encoder(x)
	pr,hh = model.decoder(zz)

	# df_z = pd.DataFrame(hh.to('cpu').detach().numpy())
	df_z = pd.DataFrame(zz.to('cpu').detach().numpy())
	df_z.columns = ["hh"+str(i)for i in df_z.columns]

	df_z["cell"] = label
	df_z["sample"] = preprocess.cellid_to_meta_single(label)
	df_z = df_z[ ["cell"]+\
			[x for x in df_z.columns if x not in["cell","sample"]]+\
			["sample"]]

	
	# data_color = range(len(df_z.columns))
	# data_color = [x / max(data_color) for x in data_color] 
	# custom_map = plt.cm.get_cmap('coolwarm') 
	# custom = custom_map(data_color)  
	# df_z.plot(kind='bar', stacked=True, color=custom,figsize=(25,10))
	# plt.ylabel("hidden state proportion", fontsize=18)
	# plt.xlabel("samples", fontsize=22)
	# plt.xticks([])
	# plt.title(title,fontsize=25)
	# plt.savefig("../output/hh_"+title+".png");plt.close()

	plt.plot(loss_values)
	plt.ylabel("loss", fontsize=18)
	plt.xlabel("epochs", fontsize=22)
	plt.title(title,fontsize=25)
	plt.savefig("../output/sc_"+title+"_loss.png");plt.close()
	df_z.to_csv("../output/sc_zz_"+title+"_data.csv",sep="\t",index=False)