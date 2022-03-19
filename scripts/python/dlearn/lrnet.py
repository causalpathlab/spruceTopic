import torch; torch.manual_seed(0)
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
import numpy as np
import pandas as pd
import annoy 
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
        indices = self.index.get_nns_by_vector(vector.tolist(),k)                                           
        return [self.labels[i] for i in indices]

class LRDataset(Dataset):
	def __init__(self,l_mat,r_mat):
		self.l_mat = l_mat
		self.r_mat = r_mat

	def __len__(self):
		return self.l_mat.shape[0]
	def __getitem__(self, idx):
		return [self.l_mat[idx], self.r_mat[idx]]


def load_data(args,nbr_size,batch_size,device):

	df_h = pd.read_csv(args.home+args.output+args.model["out"]+args.model["mfile"]+"etm_hh_data.tsv",sep="\t")

	fname = args.home+args.input+args.raw_data
	l_fname = fname.replace(".pkl",".ligands.pkl")
	r_fname = fname.replace(".pkl",".receptors.pkl")
	lr_fname = fname.replace(".pkl",".ligands_receptors_mat_815_780.pkl")

	df_l = pd.read_pickle(l_fname)
	df_r = pd.read_pickle(r_fname)
	
	df_lr = pd.read_pickle(lr_fname)
	df_lr = df_lr.loc[df_l.columns[1:],df_r.columns[1:]]
	A_lr = torch.tensor(df_lr.values.astype(np.int16))

	model_ann = ApproxNN(df_h.iloc[:,1:].to_numpy(), df_h.iloc[:,0].values)
	model_ann.build()

	cell = df_h.iloc[0,0]
	neighbours = model_ann.query(df_h.iloc[0,1:],k=nbr_size)
	l_mat,r_mat = [],[]
	for n in neighbours[1:]:
		c_m_l = torch.tensor(df_l[df_l["index"]==cell].values[0][1:].astype(np.int16)).reshape(1,df_l.shape[1]-1)
		c_n_l = torch.tensor(df_l[df_l["index"]==n].values[0][1:].astype(np.int16)).reshape(1,df_l.shape[1]-1)
		c_m_r = torch.tensor(df_r[df_r["index"]==cell].values[0][1:].astype(np.int16)).reshape(1,df_r.shape[1]-1)
		c_n_r = torch.tensor(df_r[df_r["index"]==n].values[0][1:].astype(np.int16)).reshape(1,df_r.shape[1]-1)

		l_mat.append(torch.mm(c_m_l,A_lr).mul(c_m_r) + 
					torch.mm(c_n_l,A_lr).mul(c_n_r))

		r_mat.append(torch.mm(c_m_r,torch.t(A_lr)).mul(c_m_l) + 
					torch.mm(c_n_r,torch.t(A_lr)).mul(c_n_l))

	l_mat = torch.stack(l_mat).to(device)
	r_mat = torch.stack(r_mat).to(device)
	return DataLoader(LRDataset(l_mat,r_mat), batch_size=batch_size, shuffle=True)

