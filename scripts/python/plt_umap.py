import os
import numpy as np
import pandas as pd
from collections import namedtuple
from gen_util.io import read_config
import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [12.50, 10.50]
plt.rcParams['figure.autolayout'] = True
import umap
import umap.plot 



def plot_umap_from_model():
	
	# spath = os.path.dirname(__file__)
	spath =  '/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/scripts/python'
	# spath =  '/home/sishirsubedi/projects/tumour_immune_interaction/scripts/python'
	args_home = spath.replace('/scripts/python','/')

	params = read_config(args_home+'/config/scmetm.yaml')
	args = namedtuple('Struct',params.keys())(*params.values())



	df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_zz_data.tsv',sep='\t',compression='gzip')

	umap_2d = umap.UMAP(n_components=2, init='random', random_state=0)
	proj_2d = umap_2d.fit(df_h.iloc[:,1:])
	# plt.scatter(proj_2d[:,0],proj_2d[:,1],s=0.001)
	
	df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_hh_data.tsv',sep='\t',compression='gzip')
	df_h.columns = [ x.replace('hh','k') for x in df_h.columns]
	topic_labels = df_h.iloc[:,1:].idxmax(axis=1)


	# from sklearn.cluster import KMeans
	# kmeans = KMeans(n_clusters=df_h.shape[1]-1, random_state=0).fit(df_h.iloc[:,1:].to_numpy())
	umap.plot.points(proj_2d,labels=topic_labels)

	plt.savefig('umap_zz_v3.png',dpi=300);plt.close()

	# ##test to get neighbors
	# embeddings = umap.UMAP(n_neighbors=15, min_dist=0.5).fit_transform(df_h.iloc[:,1:])
	# knn = umap.umap_.nearest_neighbors(embeddings, 
	#         n_neighbors=15, metric='euclidean', 
	#         metric_kwds={}, angular=True, 
	#         random_state=np.random.RandomState(42))

	# import pickle

	# f_name = 'umap_model.pkl'
	# pickle.dump(umap_2d, open(f_name, 'wb'))

def label_pcs():

	# spath = os.path.dirname(__file__)
	spath =  '/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/scripts/python'
	# spath =  '/home/sishirsubedi/projects/tumour_immune_interaction/scripts/python'
	args_home = spath.replace('/scripts/python','/')

	params = read_config(args_home+'/config/scmetm.yaml')
	args = namedtuple('Struct',params.keys())(*params.values())


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