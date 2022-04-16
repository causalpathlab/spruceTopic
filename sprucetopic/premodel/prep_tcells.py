import os
import pandas as pd
import datatable as dt
from scipy import sparse 
from premodel import processing as prc 
import numpy as np
import logging
logger = logging.getLogger(__name__)

def tcell_read_data_to_sparse(args,maxcellcount_pertissue=1e5,mingene_totalcount=10):
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
		
		df = prc.filter_minimal(logger,args,df,'dataset', mingene_totalcount)
		
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
	prc.tcell_lr_exp_data(args,df_combine)

	## drop lr for nbr model
	logger.info('before removing lr..'+str(df_combine.shape[0])+'x'+str(df_combine.shape[1]))
	ligand_receptors = prc.get_ligand_receptor_genes(args)
	drop_lr_columns = [x for x in df_combine.columns if x in ligand_receptors]
	df_combine = df_combine.drop(drop_lr_columns,axis=1)
	logger.info('after removing lr..'+str(df_combine.shape[0])+'x'+str(df_combine.shape[1]))

	## generate npz files
	logger.info('processing--creating coo sparse matrix file')

	label = df_combine.iloc[:,0].to_numpy()
	np.savez_compressed(sample_procpath+args.nbr_model['sparse_label'], idx=label,allow_pickle=True)

	df_combine = df_combine.iloc[:,1:]
	immune_signatures = prc.tcell_signature_genes(args_home,args)
	immune_index = [(x,y) for x,y in enumerate(df_combine.columns) if y in immune_signatures]
	np.savez_compressed(sample_procpath + args.nbr_model['sparse_immune_label'], idx=immune_index, allow_pickle=True)
	pd.DataFrame(df_combine.columns).to_pickle(sample_procpath+args.raw_data_no_lr_genes)

	S = sparse.coo_matrix(df_combine.to_numpy())
	idx, idy, val = sparse.find(S)
	d = df_combine.shape
	np.savez_compressed(sample_procpath + args.nbr_model['sparse_data'], idx=idx,idy=idy,val=val,shape=d,allow_pickle=True)

	logger.info('Data pre-processing--COMPLETED !!')

def tcell_cellid_to_meta_single(cell_ids,args):
	args_home = os.environ['args_home'] 
	meta_path = args_home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep='\t')
	df_meta = df_meta[['cellID','meta.cluster']]

	df = pd.DataFrame(cell_ids,columns=['cellID'])

	dfjoin = pd.merge(df,df_meta,on='cellID',how='left')

	return dfjoin['meta.cluster'].values


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
	