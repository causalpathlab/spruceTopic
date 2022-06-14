from torch.utils.data import DataLoader
from torch.utils.data import Dataset
import torch
from scipy import sparse
import numpy as np
import logging
logger = logging.getLogger(__name__)


class SparseData():
	def __init__(self,indptr,indices,vals,shape,label):
		self.indptr = indptr
		self.indices = indices
		self.vals = vals
		self.shape = shape
		self.label = label

class SparseDataset(Dataset):
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

def load_data(sparse_arr,arr_id,device,bath_size):

	logger.info('loading sparse data...\n'+
		sparse_arr +'\n'+ arr_id )

	npzarrs = np.load(sparse_arr,allow_pickle=True)
	s = sparse.csr_matrix( (npzarrs['val'].astype(np.int32), (npzarrs['idx'], npzarrs['idy']) ),shape=npzarrs['shape'] )

	device = torch.device(device)
	indptr = torch.tensor(s.indptr.astype(np.int32), dtype=torch.int32, device=device)
	indices = torch.tensor(s.indices.astype(np.int32), dtype=torch.int32, device=device)
	vals = torch.tensor(s.data.astype(np.int32), dtype=torch.int32, device=device)
	shape = tuple(npzarrs['shape'])

	ids = np.load(arr_id,allow_pickle=True)['idx']

	spdata = SparseData(indptr,indices,vals,shape,ids)

	return DataLoader(SparseDataset(spdata,device), batch_size=bath_size, shuffle=True)
