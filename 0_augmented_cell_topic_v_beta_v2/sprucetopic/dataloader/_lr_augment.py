from torch.utils.data import DataLoader
from torch.utils.data import Dataset
import pytorch_lightning as pl
import torch
import pandas as pd
import numpy as np

class LRDataset(Dataset):
	def __init__(self,lmat,rmat,Alr,nbrmat) :
		self.lmat = lmat
		self.rmat = rmat
		self.lrmat = Alr
		self.nbrmat = nbrmat

	def __len__(self):
		return len(self.nbrmat)

	def __getitem__(self, idx):

		cm_l = self.lmat[idx].unsqueeze(0)
		cm_r = self.rmat[idx].unsqueeze(0)

		cm_lr = torch.mm(cm_l,self.lrmat).mul(cm_r)
		cm_rl = torch.mm(cm_r,torch.t(self.lrmat)).mul(cm_l)

		nbr_idxs = self.nbrmat[idx]
		
		cn_l = self.lmat[nbr_idxs]
		cn_r = self.rmat[nbr_idxs]

		return 	cm_lr + torch.mm(cn_l,self.lrmat).mul(cn_r) , cm_rl + torch.mm(cn_r,torch.t(self.lrmat)).mul(cn_l)

class load_data(pl.LightningDataModule):

	def __init__(self,f_latent_h,f_l,f_r,f_lr,f_neighbour, batch_size,device):
		super().__init__()
		self.f_latent_h = f_latent_h
		self.f_l = f_l
		self.f_r = f_r
		self.f_lr = f_lr
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

		dflatent = pd.merge(df_l['index'],df_h,how='left',left_on='index',right_on='cell')
		dflatent['cell'] = dflatent.iloc[:,2:].idxmax(axis=1)
		dflatent = dflatent.rename(columns={'cell':'topic'})

		lmat = torch.tensor(df_l.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(self.device)
		rmat = torch.tensor(df_r.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(self.device)

		df_lr = pd.read_pickle(self.f_lr)
		df_lr = df_lr.loc[df_l.columns[1:],df_r.columns[1:]]
		Alr = torch.tensor(df_lr.values.astype(np.float32),requires_grad=False).to(self.device)

		return DataLoader(LRDataset(lmat,rmat,Alr,nbrmat), batch_size=self.batch_size, shuffle=True)

