import os
import pandas as pd
from scipy import sparse
from premodel import processing as prc 
import numpy as np
import logging
logger = logging.getLogger(__name__)

def pbmc_scanpy_processing(args):

	import scanpy as sc
	
	args_home = os.environ['args_home'] 
	adata = sc.read_10x_mtx(
			args_home+args.input+args.raw_data_scpy,  # the directory with the `.mtx` file
			var_names='gene_symbols', # use gene symbols for the variable names 
			cache=True)    

	df = adata.to_df()
	df.columns = adata.var_names
	df = df.reset_index()

	df = prc.filter_minimal(df,3)

	logger.info('fill nans as zero...')
	df.values[df.isna()] = 0

	logger.info('save raw data...')
	df.to_pickle(args_home+args.input+args.raw_data)
	pd.DataFrame(df.columns[1:]).to_pickle(args_home+args.input+args.raw_data_genes)

	logger.info('save lr...')
	prc.lr_exp_data(args,df)

	## drop lr for nbr model
	logger.info('before removing lr..'+str(df.shape[0])+'x'+str(df.shape[1]))
	ligand_receptors = prc.get_ligand_receptor_genes(args)
	drop_lr_columns = [x for x in df.columns if x in ligand_receptors]
	df = df.drop(drop_lr_columns,axis=1)
	logger.info('after removing lr..'+str(df.shape[0])+'x'+str(df.shape[1]))

	## generate npz files
	logger.info('processing--creating coo sparse matrix file')

	label = df.iloc[:,0].to_numpy()
	np.savez_compressed(args_home+args.input+args.nbr_model['sparse_label'], idx=label,allow_pickle=True)

	df = df.iloc[:,1:]
	immune_signatures = prc.get_immune_genes(args)
	immune_index = [(x,y) for x,y in enumerate(df.columns) if y in immune_signatures]
	np.savez_compressed(args_home + args.input+args.nbr_model['sparse_immune_label'], idx=immune_index, allow_pickle=True)
	pd.DataFrame(df.columns).to_pickle(args_home+args.raw_data_no_lr_genes)

	S = sparse.coo_matrix(df.to_numpy())
	idx, idy, val = sparse.find(S)
	d = df.shape
	np.savez_compressed(args_home + args.input + args.nbr_model['sparse_data'], idx=idx,idy=idy,val=val,shape=d,allow_pickle=True)

	logger.info('Data pre-processing--COMPLETED !!')

