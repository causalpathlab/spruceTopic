import os
import pandas as pd
import matplotlib.pylab as plt
import sprucetopic.premodel.processing as prep
import numpy as np
import scanpy as sc

plt.rcParams['figure.figsize'] = [12.50, 10.50]
plt.rcParams['figure.autolayout'] = True

def scanpy_filter(args):
	
	args_home = os.environ['args_home'] 
	
	df = pd.read_pickle(args_home+args.input+args.raw_data)

	## convert df to anndata object
	obs = pd.DataFrame(df.iloc[:,0])
	var = pd.DataFrame(index=df.columns[1:])
	X = df.iloc[:,1:].to_numpy()
	adata = sc.AnnData(X, obs=obs, var=var, dtype='int32')

	## filter scanpy default settings
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

	## scanpy pca check
	sc.tl.pca(adata, svd_solver='arpack')
	sc.pl.pca(adata, color='CST3')
	plt.savefig('../../output/scanpy_raw_pipeline_pca.png');plt.close()
	
	sc.pl.pca_variance_ratio(adata, n_pcs=50,log=True)
	plt.savefig('../../output/scanpy_raw_pipeline_pca_var.png');plt.close()

	###scanpy results
	pd.DataFrame(adata.var.highly_variable.index).to_csv(
		'../../output/scanpy_raw_pipeline_1k_hvariable_genes.tsv',
		sep='\t',index=False, compression='gzip'
	)

	df_scanpy = adata.to_df()
	df_scanpy.columns = adata.var_names
	df_scanpy['index'] = adata.obs['index']
	df_scanpy = df_scanpy[['index']+[x for x in df_scanpy.columns[:-1]]]

	df_scanpy.to_pickle('../../output/scanpy_raw_pipeline_166k_1k_hvariable_genes.pkl')

	df_scanpy_pcs = pd.DataFrame(adata.obsm['X_pca'])
	df_scanpy_pcs['index'] = adata.obs['index'].values
	df_scanpy_pcs = df_scanpy_pcs[['index']+[int(x) for x in df_scanpy_pcs.columns[:-1]]]
	df_scanpy_pcs.to_pickle('../../output/scanpy_raw_pipeline_166k_50PCS_hvariable_genes.pkl')

def neighbour_analysis(adata):
	sc.pp.neighbors(adata)
	sc.tl.umap(adata)
	sc.pl.umap(adata, color=['CST3'])
	plt.savefig('../output/scanpy_etm_tcell_all_umap.png');plt.close()

	sc.tl.leiden(adata)
	sc.pl.umap(adata, color=['leiden'])
	plt.savefig('../output/scanpy_tcell_all_umap_leiden.png');plt.close()

	pclust = prep.cellid_to_meta_single(adata.obs['index'].values)

	adata.obs['pclust']=pclust
	sc.pl.umap(adata, size=2.0, color=['pclust'])
	plt.savefig('../output/scanpy_tcell_all_umap_pclust.png');plt.close()

	marker = ['CD4.c20.Treg.TNFRSF9',\
		'CD4.c06.Tm.ANXA1',\
		'CD4.c01.Tn.TCF7',\
		'CD8.c02.Tm.IL7R',\
		'CD8.c05.Tem.CXCR5',\
		'CD8.c07.Temra.CX3CR1']

	pclust_mod = [ x if x in marker else 'others' for x in pclust]
	adata.obs['pclust_mod']=pclust_mod
	clrs = ['dodgerblue','purple','green','red','orange','brown','grey']
	sc.pl.umap(adata, size=2.0, color=['pclust_mod'],palette=clrs)
	plt.savefig('../output/scanpy_etm_tcell_all_umap_pclust_selected.png');plt.close()


	pclust_mod2 = [ 'CD4' if 'CD4.' in x else 'CD8' for x in pclust]
	adata.obs['pclust_mod2']=pclust_mod2
	clrs = ['dodgerblue','orange',]
	sc.pl.umap(adata, size=2.0, color=['pclust_mod2'],palette=clrs)
	plt.savefig('../output/scanpy_tcell_all_umap_pclust_selected_CD4_8.png');plt.close()
	# component_check.run_component_analysis(
	#     ,
	#     ,
	#     'tsne',data.split('.')[0])
	# component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,'pca',data.split('.')[0])
 
	## get df back from anndata object
	df_scanpy = adata.to_df()
	df_scanpy.columns = adata.var_names
	df_scanpy['cell'] = adata.obs['index']
	df_scanpy['leiden'] = adata.obs['leiden']
	df_scanpy['pclust'] = adata.obs['plust']
	df_scanpy['pclust_mod'] = adata.obs['plust_mod']

	import umap
	import seaborn as sns
	reducer = umap.UMAP()
	embedding = reducer.fit_transform(df_scanpy.iloc[:,:-2])
	clrs = [sns.color_palette()[int(x)] for x in df_scanpy['leiden'].unique()]
	plt.scatter(
	embedding[:, 0],
	embedding[:, 1],s=0.1)
	plt.gca().set_aspect('equal', 'datalim')
	plt.title('UMAP projection', fontsize=24)
	plt.savefig('../output/scanpy_tcell_all_umap_from_library.png');plt.close()


def label_pcs(args):
	import umap
	import umap.plot 

	args_home = os.environ['args_home'] 
	df = pd.read_pickle('../../output/scanpy_raw_pipeline_166k_50PCS_hvariable_genes.pkl')

	meta_path = args_home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep='\t')
	labels = pd.merge(df[['index']],df_meta,right_on='cellID',left_on='index',how='left')['cancerType'].values

	umap_2d = umap.UMAP(n_components=2, init='random', random_state=0)
	proj_2d = umap_2d.fit(df.iloc[:,1:])
	umap.plot.points(proj_2d,labels=labels)
	plt.savefig('scanpy_50pcs_celltype_umap.png',dpi=300);plt.close()


	labels = pd.merge(df[['index']],df_meta,right_on='cellID',left_on='index',how='left')['meta.cluster'].values

	cdlabels = np.array([x.split('.')[0] for x in labels])
	umap.plot.points(proj_2d,labels=cdlabels)
	plt.savefig('scanpy_50pcs_cellcluster_cd_umap.png',dpi=300);plt.close()

	marker = ['CD4.c20.Treg.TNFRSF9',\
		'CD4.c06.Tm.ANXA1',\
		'CD4.c01.Tn.TCF7',\
		'CD8.c02.Tm.IL7R',\
		'CD8.c05.Tem.CXCR5',\
		'CD8.c07.Temra.CX3CR1']

	cdlabels_top = np.array([ x if x in marker else 'others' for x in labels])
	umap.plot.points(proj_2d,labels=cdlabels_top)
	plt.savefig('scanpy_50pcs_cellcluster_top_umap.png',dpi=300);plt.close()


	df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_hh_data.tsv',sep='\t',compression='gzip')
	df_h.columns = [ x.replace('hh','k') for x in df_h.columns]

	topic_labels = df_h.iloc[:,1:].idxmax(axis=1)
	df_h['label'] = topic_labels
	df_h = df_h[['cell','label']]

	etm_labels = pd.merge(df[['index']],df_h,right_on='cell',left_on='index',how='left')['label'].values
	umap.plot.points(proj_2d,labels=etm_labels)
	plt.savefig('scanpy_50pcs_cellcluster_etm_umap.png',dpi=300);plt.close()


def pbmc_scanpy_processing(args):
	
	args_home = os.environ['args_home'] 
	adata = sc.read_10x_mtx(
			args_home+args.input+args.raw_data_scpy,  # the directory with the `.mtx` file
			var_names='gene_symbols', # use gene symbols for the variable names 
			cache=True)    

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
	sc.pl.pca(adata, color='CST3')
	plt.savefig('pbmc_scanpy_raw_pipeline_pca.png');plt.close()
	
	sc.pl.pca_variance_ratio(adata, n_pcs=50,log=True)
	plt.savefig('pbmc_scanpy_raw_pipeline_pca_var.png');plt.close()

	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
	sc.tl.umap(adata)
	sc.tl.leiden(adata)
	sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
	plt.savefig('pbmc_scanpy_raw_pipeline_umap.png');plt.close()



	

