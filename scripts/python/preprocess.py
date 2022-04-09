import os
import pandas as pd
import matplotlib.pylab as plt
from gen_util.io import read_config
import datatable as dt
from scipy import sparse 
import numpy as np
import logging
logger = logging.getLogger(__name__)
plt.rcParams['figure.figsize'] = [12.50, 10.50]
plt.rcParams['figure.autolayout'] = True


def select_protein_coding_genes(df):
	params = read_config()
	pcg = params['home']+params['database']+params['protein_coding_genes']
	df_pcg = pd.read_csv(pcg,header=None,names=['gene'])
	drop_columns = [x for x in df.columns[1:-1] if x not in df_pcg['gene'].values]
	df = df.drop(drop_columns,axis=1)
	return df

def get_ligand_receptor_genes(args):
	df = pd.read_csv( os.environ['args_home']+args.database+args.lr_db,sep='\t')
	return list(df.receptor_gene_symbol.unique())+list(df.ligand_gene_symbol.unique())

def filter_minimal(logger,args,df,mode, cutoff):
	
	logger.info('initial data size : '+str(df.shape[0])+'x'+str(df.shape[1]))

	drop_columns =[]

	ligand_receptors = get_ligand_receptor_genes(args)
	drop_columns =[]

	if mode == 'tissue':
		# eliminate gene if the total number of counts is < cutoff per tissue type.
		for tissue in df['sample'].unique():
			drop_columns_sample = [ col for col,val  in df[df['sample']==tissue].iloc[:,1:-1].sum(axis=0).iteritems() if val < cutoff ]
			for c in drop_columns_sample:
				if c not in drop_columns:
					drop_columns.append(c)
	
	elif mode == 'dataset':
		drop_columns = [ col for col,val  in df.iloc[:,1:].sum(axis=0).iteritems() if val < cutoff ]

	logger.info('genes to filter based on mincout cutoff - '+ str(len(drop_columns)))

	for mt_g in [x for x in df.columns if 'MT-' in x]:
		drop_columns.append(mt_g)

	logger.info('adding mitochondrial genes - '+ str(len(drop_columns)))

	for spk_g in [x for x in df.columns if 'ERCC' in x]:
		drop_columns.append(mt_g)

	logger.info('adding spikes - '+ str(len(drop_columns)))

	df = df.drop(drop_columns,axis=1)
	logger.info('after all minimal filtering..'+str(df.shape[0])+'x'+str(df.shape[1]))

	return df

def read_sc_data_sparse(args,maxcellcount_pertissue=1e5,mingene_totalcount=10):
	args_home = os.environ['args_home'] 
	df_combine = pd.DataFrame()
	sample_path = args_home+args.input+args.nbr_model['in']
	sample_procpath = args_home+args.input

	for sample in args.samples:
		
		logger.info('processing--'+sample)
		
		dtab = dt.fread(sample_path+sample,verbose=False)
		df = dtab.to_pandas() 
		df = df.T
		df = df.rename(columns=df.iloc[0])
		df = df.iloc[1:].reset_index()
		
		df = filter_minimal(logger,args,df,'dataset', mingene_totalcount)
		
		if df.shape[0]>maxcellcount_pertissue:
			df = df.sample(n=maxcellcount_pertissue)
		if sample == args.sample_first:
			df_combine = df
			logger.info('current: '+str(df.shape[0])+'x'+str(df.shape[1]))
			logger.info('combine: '+str(df_combine.shape[0])+'x'+str(df_combine.shape[1]))
		else:
			logger.info('concatinating...')
			df_combine = pd.concat([df_combine, df], axis=0, ignore_index=True)
			logger.info('current: '+str(df.shape[0])+'x'+str(df.shape[1]))
			logger.info('combine: '+str(df_combine.shape[0])+'x'+str(df_combine.shape[1]))
	
	del df
	logger.info('fill nans as zero...')
	df_combine.values[df_combine.isna()] = 0

	logger.info('save raw data...')
	df_combine.to_pickle(sample_procpath+args.raw_data)
	pd.DataFrame(df_combine.columns[1:]).to_pickle(sample_procpath+args.raw_data_genes)

	logger.info('save lr...')
	create_lr_exp_data(args,df_combine)

	## drop lr for nbr model
	logger.info('before removing lr..'+str(df_combine.shape[0])+'x'+str(df_combine.shape[1]))
	ligand_receptors = get_ligand_receptor_genes(args)
	drop_lr_columns = [x for x in df_combine.columns if x in ligand_receptors]
	df_combine = df_combine.drop(drop_lr_columns,axis=1)
	logger.info('after removing lr..'+str(df_combine.shape[0])+'x'+str(df_combine.shape[1]))

	## generate npz files
	logger.info('processing--creating coo sparse matrix file')

	label = df_combine.iloc[:,0].to_numpy()
	np.savez_compressed(sample_procpath+args.nbr_model['sparse_label'], idx=label,allow_pickle=True)

	df_combine = df_combine.iloc[:,1:]
	immune_signatures = get_tcell_signature_genes(args_home,args)
	immune_index = [(x,y) for x,y in enumerate(df_combine.columns) if y in immune_signatures]
	np.savez_compressed(sample_procpath + args.nbr_model['sparse_immune_label'], idx=immune_index, allow_pickle=True)

	S = sparse.coo_matrix(df_combine.to_numpy())
	idx, idy, val = sparse.find(S)
	d = df_combine.shape
	np.savez_compressed(sample_procpath + args.nbr_model['sparse_data'], idx=idx,idy=idy,val=val,shape=d,allow_pickle=True)

	logger.info('Data pre-processing--COMPLETED !!')


def select_cells():
	'''
	Create a test data for algorithm development from Zheng et al paper (GSE156728) 
	selected data - tissue type, cell type, cluster.
	'''
	
	celltype_cluster = ['CD4_c18','CD4_c13']
	tissue_type = ['BC', 'BCL', 'ESCA', 'MM', 'PACA', 'RC', 'THCA', 'UCEC']
	
	params = read_config()
	meta_path = params['home']+params['database']
	
	df = pd.read_csv(meta_path+params['metadata'],sep='\t')
	
	df['meta.cluster.sf'] =[x.split('.')[0]+'_'+x.split('.')[1] for x in df['meta.cluster']]
	df = df[(df['meta.cluster.sf'].isin(celltype_cluster))]
	df = df[df['cancerType'].isin(tissue_type)]
	df['cluster'] = [x+'_'+y for x,y in zip(df['cancerType'],df['meta.cluster.sf'])]
	df = df[['cellID','cluster']]
	df.to_csv('../output/all_tissue_cd4_c13_c18_cell_barcodes.csv',sep='\t',index=False)

def cellid_to_meta(cell_pairs,args):
	meta_path = args.home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep='\t')

	cell_pairs_meta = []
	for cp in cell_pairs:
		types = df_meta.loc[df_meta['cellID'].isin(cp),['meta.cluster']].values
		cell_pairs_meta.append(types.flatten()[0]+'/'+types.flatten()[1])

	return cell_pairs_meta

def cellid_to_meta_single(cell_ids,args):
	args_home = os.environ['args_home'] 
	meta_path = args_home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep='\t')
	df_meta = df_meta[['cellID','meta.cluster']]

	df = pd.DataFrame(cell_ids,columns=['cellID'])

	dfjoin = pd.merge(df,df_meta,on='cellID',how='left')

	return dfjoin['meta.cluster'].values

def get_immune_genes():
	params = read_config()
	meta_path = params['home']+params['database']
	df_meta = pd.read_csv(meta_path+params['immune_signature_genes'],sep='\t')
	return df_meta['signature_genes'].values

def get_tcell_signature_genes(args_home,args):
	df_meta = pd.read_csv(args_home+args.database+'pancancer_tcell_signature_gene_all_clusters.tsv')
	return df_meta['signature_genes'].values

def create_lr_exp_data(args,df_exp):
	args_home = os.environ['args_home']
	df = pd.read_csv( args_home + args.database+args.lr_db,sep='\t')
	receptors = list(df.receptor_gene_symbol.unique())
	ligands = list(df.ligand_gene_symbol.unique())

	l_fname = args_home+args.input+args.raw_l_data
	r_fname = args_home+args.input+args.raw_r_data

	df_exp[['index']+[ x for x in ligands if x in df_exp.columns ]].to_pickle(l_fname)
	df_exp[['index']+[ x for x in receptors if x in df_exp.columns ]].to_pickle(r_fname)

def create_lr_mat(args):
	df = pd.read_csv( os.environ['args_home']+args.database+args.lr_db,sep='\t')
	dflrmat = df.groupby(['ligand_gene_symbol','receptor_gene_symbol']).agg(['count'])['lr_pair']
	dflrmat = dflrmat.unstack(fill_value=0)
	dflrmat.columns = dflrmat.columns.droplevel(0)
	fname = os.environ['args_home']+args.input+args.raw_lr_data
	dflrmat.to_pickle(fname)

def tcell_marker_genes(args_home,args,tag,mode):

	db = args_home+args.database+args.tcell_signature_genes

	if mode =='sep':
		df = pd.read_excel(db,sheet_name=tag,header=1,usecols='A,B,C' )
		# df_cd8 = pd.read_excel(db,sheet_name='CD8',header=1,usecols='A,B,C' )
		# df = df_cd4.append(df_cd8,ignore_index=True)
		df['clust'] = [x.split('(')[0] for x in df['cluster.name']]
		df = df[['geneSymbol','clust']]
		df_gene_clust = df.groupby(['geneSymbol'], as_index = False).agg('/'.join)
		df_gene_clust['clust_id'] = pd.factorize(df_gene_clust['clust'])[0]
		df_gene_clust.to_csv(args_home+args.database+'tcell_signature_gene_clusters_'+tag+'.tsv',sep='\t',index=False)
	
	elif mode =='combine':

		df_cd4 = pd.read_excel(db,sheet_name='CD4',header=1,usecols='A,B,C')
		df_cd8 = pd.read_excel(db,sheet_name='CD8',header=1,usecols='A,B,C')

		df_tcell_markers = pd.DataFrame(pd.Series(np.concatenate([df_cd4['geneSymbol'].values,df_cd8['geneSymbol'].values])).unique())
		df_tcell_markers = df_tcell_markers.sort_values(0)
		df_tcell_markers.to_csv(args_home+args.database+'pancancer_tcell_signature_gene_all_clusters.tsv',sep='\t',index=False)