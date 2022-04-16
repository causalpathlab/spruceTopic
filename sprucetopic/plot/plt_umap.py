import os
import numpy as np
import pandas as pd
from collections import namedtuple
from _utils.io import read_config
import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [12.50, 10.50]
plt.rcParams['figure.autolayout'] = True
import umap
import umap.plot 

def plot_umap_from_model(args):

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
	df_z = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_zz_data.tsv',sep='\t',compression='gzip')

	df_meta = pd.read_pickle(args_home+args.scanpy_output+'pbmc_cluster.pkl')
	dfjoin = pd.merge(df_z,df_meta,right_on='index',left_on='cell',how='left')
	dfjoin = dfjoin[dfjoin['leiden'].notna()]
	
	umap_2d = umap.UMAP(n_components=2, init='random', random_state=0,min_dist=0.1)
	proj_2d = umap_2d.fit(dfjoin.iloc[:,1:-2])
	celltype_labels = np.array([ cell_type[int(x)] for x in dfjoin['leiden'].values ])
	print(celltype_labels)
	umap.plot.points(proj_2d,labels=celltype_labels,color_key_cmap='Paired')
	plt.savefig(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'umap_zz_scanpy_cluster_label.png',dpi=300);plt.close()

	
	df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_hh_data.tsv',sep='\t',compression='gzip')
	df_h.columns = [ x.replace('hh','k') for x in df_h.columns]
	df_h['topic'] = df_h.iloc[:,1:].idxmax(axis=1)

	topic_labels = pd.merge(dfjoin[['cell']],df_h[['cell','topic']],on='cell',how='left')['topic'].values

	umap.plot.points(proj_2d,labels=topic_labels,color_key_cmap='tab20b')
	plt.savefig(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'umap_zz_topic_label.png',dpi=300);plt.close()




def label_pcs(args):
	
	args_home = os.environ['args_home'] 

	df = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_zz_data.tsv',sep='\t',compression='gzip')

	meta_path = args_home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep='\t')
	labels = pd.merge(df[['cell']],df_meta,right_on='cellID',left_on='cell',how='left')['cancerType'].values

	umap_2d = umap.UMAP(n_components=2, init='random', random_state=0)
	proj_2d = umap_2d.fit(df.iloc[:,1:])

	umap.plot.points(proj_2d,labels=labels)
	plt.savefig('etm_50pcs_celltype_umap.png', dpi=300, format='png')
	plt.close()


	labels = pd.merge(df[['cell']],df_meta,right_on='cellID',left_on='cell',how='left')['meta.cluster'].values

	cdlabels = np.array([x.split('.')[0] for x in labels])
	umap.plot.points(proj_2d,labels=cdlabels)
	plt.savefig('etm_50pcs_cellcluster_cd_umap.png',dpi=300);plt.close()

	marker = ['CD4.c20.Treg.TNFRSF9',\
		'CD4.c06.Tm.ANXA1',\
		'CD4.c01.Tn.TCF7',\
		'CD8.c02.Tm.IL7R',\
		'CD8.c05.Tem.CXCR5',\
		'CD8.c07.Temra.CX3CR1']

	cdlabels_top = np.array([ x if x in marker else 'others' for x in labels])
	umap.plot.points(proj_2d,labels=cdlabels_top)
	plt.savefig('etm_50pcs_cellcluster_top_umap.png',dpi=300);plt.close()