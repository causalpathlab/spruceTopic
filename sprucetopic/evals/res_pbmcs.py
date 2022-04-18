import os
import pandas as pd
import numpy as np

def pbmc_sample_cells_with_latent(args,cell_n=50):

	cell_type={
		0:'CD4T-Naive',
		1:'CD4T-Memory',
		4:'CD8T',
		6:'NK',
		3:'B',
		2:'CD14',
		5:'FCGR3A',
		7:'Dendritic',
		8:'Megakaryocytes'
	}

	args_home = os.environ['args_home']

	df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_hh_data.tsv',sep='\t',compression='gzip')

	meta_path = args_home+args.scanpy_output+'pbmc_cluster.pkl'
	df_meta = pd.read_pickle(meta_path)
	dfjoin = pd.merge(df_h,df_meta,right_on='index',left_on='cell',how='left')
	dfjoin = dfjoin[dfjoin['leiden'].notna()]


	# df_h_sample = dfjoin.groupby('leiden').sample(cell_n, random_state=1)
	df_h_sample = dfjoin.drop(columns='index').rename(columns={'leiden':'cluster'})
	df_h_sample['cluster'] = [cell_type[int(x)] for x in df_h_sample['cluster']]
	
	df_h_sample.columns = [ x.replace('hh','k') for x in df_h_sample.columns]
	df_h_sample.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'hh_cell_topic_sample.tsv',sep='\t',index=False)



def topic_celltype_marker_genes(args):

	args_home = os.environ['args_home']

	df_genes=pd.read_pickle(args_home+args.input+args.raw_data_no_lr_genes)

	immuneindex_file = args_home+args.input+args.nbr_model['sparse_immune_label']
	immunedat = np.load(immuneindex_file,allow_pickle=True)
	immuneindx = np.array([x[0] for x in immunedat['idx']]).astype(np.int32)
	otherindx = np.array([x for x in range(df_genes.shape[0]) if x not in immuneindx]).astype(np.int32)

	df_beta1 = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_beta1_data.tsv',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_beta2_data.tsv',sep='\t',compression='gzip')

	df_beta1.columns = df_genes.iloc[immuneindx][0].values
	df_beta2.columns = df_genes.iloc[otherindx][0].values

	marker_genes = ['IL7R', 'S100A4','CCR7','CD79A', 'MS4A1', 'CD8A', 'CD8B', 
	               'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
	
	# In [43]: [x for x in ligand_receptors if x in marker_genes]
	# Out[43]: ['IL7R', 'FCER1A', 'CD79A', 'KLRB1', 'LYZ', 'PPBP', 'CD14', 'S100A8']
	
	df_beta = pd.concat([df_beta1,df_beta2],axis=1)

	marker_genes = [x for x in marker_genes if x in df_beta.columns]

	df_beta = df_beta[marker_genes]

	df_beta.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'marker_genes_topic.tsv',sep='\t',index=False)

def topic_celltype_marker_genes_lrnet(args):

	args_home = os.environ['args_home']

	df_l_genes=pd.read_pickle(args_home+args.input+args.raw_l_data_genes)
	df_r_genes=pd.read_pickle(args_home+args.input+args.raw_r_data_genes)


	df_beta1 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta1_data.tsv',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta2_data.tsv',sep='\t',compression='gzip')

	df_beta1.columns = df_r_genes[0].values
	df_beta2.columns = df_l_genes[0].values

	marker_genes = ['IL7R', 'S100A4','CCR7','CD79A', 'MS4A1', 'CD8A', 'CD8B', 
	               'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
	
	# In [43]: [x for x in ligand_receptors if x in marker_genes]
	# Out[43]: ['IL7R', 'FCER1A', 'CD79A', 'KLRB1', 'LYZ', 'PPBP', 'CD14', 'S100A8']
	
	df_beta = pd.concat([df_beta1,df_beta2],axis=1)

	marker_genes = [x for x in marker_genes if x in df_beta.columns]

	df_beta = df_beta[marker_genes]

	df_beta.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'marker_genes_topic.tsv',sep='\t',index=False)
