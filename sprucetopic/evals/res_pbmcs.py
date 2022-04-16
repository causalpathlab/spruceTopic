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



