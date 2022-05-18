import os
import pandas as pd
import datatable as dt
import numpy as np
import logging
logger = logging.getLogger(__name__)


def select_protein_coding_genes(args,df):
	df_pcg = pd.read_csv(os.environ['args_home']+args.database+args.protein_coding_genes,header=None,names=['gene'])
	drop_columns = [x for x in df.columns[1:-1] if x not in df_pcg['gene'].values]
	df = df.drop(drop_columns,axis=1)
	return df

def get_ligand_receptor_genes(args):
	df = pd.read_csv( os.environ['args_home']+args.database+args.lr_db,sep='\t')
	return list(df.receptor_gene_symbol.unique())+list(df.ligand_gene_symbol.unique())

def filter_minimal(df,cutoff):
	
	logger.info('initial data size : '+str(df.shape[0])+'x'+str(df.shape[1]))

	drop_columns = [ col for col,val  in df.iloc[:,1:].sum(axis=0).iteritems() if val < cutoff ]

	logger.info('genes to filter based on mincout cutoff - '+ str(len(drop_columns)))

	for mt_g in [x for x in df.columns if 'MT-' in x]:
		drop_columns.append(mt_g)

	logger.info('adding mitochondrial genes - '+ str(len(drop_columns)))

	for spk_g in [x for x in df.columns if 'ERCC' in x]:
		drop_columns.append(spk_g)

	logger.info('adding spikes - '+ str(len(drop_columns)))

	df = df.drop(drop_columns,axis=1)

	logger.info('after all minimal filtering..'+str(df.shape[0])+'x'+str(df.shape[1]))

	return df

def get_immune_genes_bcancer(args):
	args_home = os.environ['args_home']
	df_meta = pd.read_csv(args_home+args.database+args.immune_signature_genes,sep="\t")

	df_meta = df_meta[
				( df_meta['species'].isin(['Mm Hs','Hs']) ) &
			    (( df_meta['organ'].isin(['Immune system']) ) |
			    ( df_meta['organ'].isin(['Mammary gland']) ))
			]
			
	return df_meta['official gene symbol'].unique()

def lr_exp_data(args,df_exp):
	args_home = os.environ['args_home']
	df = pd.read_csv( args_home + args.database+args.lr_db,sep='\t')
	receptors = list(df.receptor_gene_symbol.unique())
	ligands = list(df.ligand_gene_symbol.unique())

	df_exp[['index']+[ x for x in ligands if x in df_exp.columns ]].to_pickle(
		args_home+args.input+args.raw_l_data)
	df_exp[['index']+[ x for x in receptors if x in df_exp.columns ]].to_pickle(
		args_home+args.input+args.raw_r_data)

	pd.DataFrame([ x for x in ligands if x in df_exp.columns ]).to_pickle(args_home+args.input+args.raw_l_data_genes)
	pd.DataFrame([ x for x in receptors if x in df_exp.columns ]).to_pickle(args_home+args.input+args.raw_r_data_genes)

def create_lr_mat(args):
	df = pd.read_csv( os.environ['args_home']+args.database+args.lr_db,sep='\t')
	dflrmat = df.groupby(['ligand_gene_symbol','receptor_gene_symbol']).agg(['count'])['lr_pair']
	dflrmat = dflrmat.unstack(fill_value=0)
	dflrmat.columns = dflrmat.columns.droplevel(0)
	fname = os.environ['args_home']+args.input+args.raw_lr_data
	dflrmat.to_pickle(fname)

def tenx_preprocessing(args):
	import scanpy as sc
	from scipy import sparse 

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
	df = filter_minimal(df,min_total_gene_count)


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

def lr_preprocessing(args):
	df = pd.read_pickle(os.environ['args_home']+args.input+args.raw_data)
	lr_exp_data(args,df)
	create_lr_mat(args)

def scanpy_processing(args):
	import scanpy as sc
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
	df_scanpy.to_pickle(args_home+ args.output+ args.sample+'scanpy_raw_pipeline_hvariable_genes.pkl')

	df_scanpy_pcs = pd.DataFrame(adata.obsm['X_pca'])
	df_scanpy_pcs['index'] = df_scanpy['index']
	df_scanpy_pcs = df_scanpy_pcs[['index']+[int(x) for x in df_scanpy_pcs.columns[:-1]]]
	df_scanpy_pcs.to_pickle(args_home+ args.output+ args.sample+'scanpy_raw_pipeline_PC_dimensions.pkl')

