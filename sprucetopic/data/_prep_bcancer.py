import os
import pandas as pd
import scanpy as sc
from scipy import sparse
from premodel import processing as prc 
import numpy as np
import logging
logger = logging.getLogger(__name__)

def preprocessing(args):

	args_home = os.environ['args_home'] 

	logger.info('reading data to scanpy...')

	adata = sc.read_10x_mtx(args_home+args.tenx_data)

	logger.info('reading data to pandas...')
	df = adata.to_df()
	df.columns = adata.var_names
	df = df.reset_index()

	# logger.info('fill nans as zero...')
	# df.values[df.isna()] = 0.0

	min_total_gene_count = 3
	df = prc.filter_minimal(df,min_total_gene_count)


	logger.info('save raw data...')
	df.to_pickle(args_home+args.input+args.raw_data)
	pd.DataFrame(df.columns[1:]).to_pickle(args_home+args.input+args.raw_data_genes)

	## generate npz files
	logger.info('processing--creating coo sparse matrix file')

	label = df.iloc[:,0].to_numpy()
	np.savez_compressed(args_home+args.input+args.nbr_model['sparse_label'], idx=label,allow_pickle=True)

	df = df.iloc[:,1:]

	S = sparse.coo_matrix(df.to_numpy())
	idx, idy, val = sparse.find(S)
	d = df.shape
	np.savez_compressed(args_home + args.input + args.nbr_model['sparse_data'], idx=idx,idy=idy,val=val,shape=d,allow_pickle=True)

	logger.info('Data pre-processing--COMPLETED !!')

def combine_metadata(args):

	import glob 

	args_home = os.environ['args_home'] 

	logger.info('reading data to scanpy...')

	path = args.meta_data

	csv_files = glob.glob(os.path.join(path, "*.csv"))

	df_combine = pd.DataFrame()

	for f in csv_files:
		df = pd.read_csv(f)
		print(f,df.shape)
		df_combine = pd.concat([df_combine, df], axis=0, ignore_index=True)

	df_combine.to_csv(args_home+ args.input+ args.sample+'_metadata.tsv.gz',sep='\t',index=False,compression='gzip')

def scanpy_processing(args):
	import matplotlib.pylab as plt
	plt.rcParams['figure.figsize'] = [12.50, 10.50]
	plt.rcParams['figure.autolayout'] = True

	
	args_home = os.environ['args_home'] 

	adata = sc.read_10x_mtx(args_home+args.tenx_data)

	sc.pp.filter_cells(adata, min_genes=200)
	sc.pp.filter_genes(adata, min_cells=3)
	adata.var['mt'] = adata.var_names.str.startswith('MT-')  
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	adata = adata[adata.obs.n_genes_by_counts < 2500, :]
	adata = adata[adata.obs.pct_counts_mt < 5, :]
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.log1p(adata)
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	adata = adata[:, adata.var.highly_variable]
	sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
	sc.pp.scale(adata, max_value=10)

	sc.tl.pca(adata, svd_solver='arpack')
	sc.pl.pca(adata)
	plt.savefig(args_home+ args.output+ args.sample+'_scanpy_raw_pipeline_pca.png');plt.close()
	
	sc.pl.pca_variance_ratio(adata, n_pcs=50,log=True)
	plt.savefig(args_home+ args.output+ args.sample+'_scanpy_raw_pipeline_pca_var.png');plt.close()

	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
	sc.tl.umap(adata)
	sc.tl.leiden(adata)
	sc.pl.umap(adata, color=['leiden'])
	plt.savefig(args_home+ args.output+ args.sample+'_scanpy_raw_pipeline_umap.png');plt.close()

	df_scanpy = adata.to_df()
	df_scanpy.columns = adata.var_names
	df_scanpy = df_scanpy.reset_index()
	df_scanpy.to_pickle(args_home+ args.output+ args.sample+'scanpy_raw_pipeline_41k_2k_hvariable_genes.pkl')

	df_scanpy_pcs = pd.DataFrame(adata.obsm['X_pca'])
	df_scanpy_pcs['index'] = df_scanpy['index']
	df_scanpy_pcs = df_scanpy_pcs[['index']+[int(x) for x in df_scanpy_pcs.columns[:-1]]]
	df_scanpy_pcs.to_pickle(args_home+ args.output+ args.sample+'scanpy_raw_pipeline_41k_50PCS_hvariable_genes.pkl')


def lr_preprocessing(args):
	
	args_home = os.environ['args_home'] 

	df = pd.read_pickle(args_home+args.input+args.raw_data)
	prc.lr_exp_data(args,df)
	prc.create_lr_mat(args)


