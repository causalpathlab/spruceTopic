from torch.utils.data import DataLoader
from torch.utils.data import Dataset
import pytorch_lightning as pl
import pandas as pd
import numpy as np
import torch


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


class load_data(pl.LightningDataModule):

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
