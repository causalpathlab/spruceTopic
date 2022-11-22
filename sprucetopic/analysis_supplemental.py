##################################################################
    # setup
##################################################################
import experiment_interaction_topic 
import pandas as pd
from pathlib import Path

server = Path.home().as_posix()
experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'
experiment_home = server+experiment_home
spr = experiment_interaction_topic.get_experiment_model(experiment_home)
print(spr.cell_topic.model_id)
print(spr.interaction_topic.model_id)


##################################################################
############ interaction topic tree
##################################################################

df = pd.read_csv(spr.interaction_topic.model_id+'_it_model_deg_analysis.csv.gz')
df_meta = pd.read_csv(spr.interaction_topic.model_id+'_alltopics_meta.csv.gz',compression='gzip')


from scipy.cluster.hierarchy import dendrogram, linkage
import sys
sys.setrecursionlimit(10000)
import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [15, 10]
plt.rcParams['figure.autolayout'] = True
import colorcet as cc
import seaborn as sns

cells = df_meta.groupby('interact_topic').sample(frac=0.1, random_state=1)['cell'].values
df_sample = df[df['cell_id'].isin(cells)]
Z = linkage(df_sample.iloc[:,1:].T.values)
dendrogram(Z)  
plt.savefig(spr.interaction_topic.model_id+'_it_hierarchy.png');plt.close()

##################################################################
############ scanpy processing of mix dataset ~ cell topic model
##################################################################

def scanpy_processing(df):
    import scanpy as sc
	import matplotlib.pylab as plt
	plt.rcParams['figure.figsize'] = [12.50, 10.50]
	plt.rcParams['figure.autolayout'] = True
	
	args_sample='/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/output/scanpy/gse_mix'
	df = df.set_index('index')
	adata = sc.AnnData(df)	

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
	plt.savefig(args_sample+'_scanpy_raw_pipeline_pca.png');plt.close()
	
	sc.pl.pca_variance_ratio(adata, n_pcs=50,log=True)
	plt.savefig(args_sample+'_scanpy_raw_pipeline_pca_var.png');plt.close()

	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
	sc.tl.umap(adata)
	sc.tl.leiden(adata)
	sc.pl.umap(adata, color=['leiden'])
	plt.savefig( args_sample+'_scanpy_raw_pipeline_umap.png');plt.close()


	df_scanpy = adata.to_df()
	df_scanpy.columns = adata.var_names
	df_scanpy = df_scanpy.reset_index()
	# df_scanpy.to_pickle(args_sample+'scanpy_raw_pipeline_123k_2.5k_hvariable_genes.pkl')


	df_scanpy_umap = pd.DataFrame(adata.obsm['X_umap'])
	df_scanpy_umap['cell'] = df_scanpy['index']
	df_umap = pd.read_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_kmeans.csv.gz')
	df_scanpy_umap = pd.merge(df_scanpy_umap,df_umap,on='cell',how='left')

	import matplotlib.pylab as plt
	plt.rcParams['figure.figsize'] = [15, 10]
	plt.rcParams['figure.autolayout'] = True
	import colorcet as cc
	import seaborn as sns
	
	cp = sns.color_palette(cc.glasbey_dark, n_colors=len(df_scanpy_umap['cluster_celltype'].unique()))
	p = sns.scatterplot(data=df_scanpy_umap, x=0, y=1, hue='cluster_celltype',s=2,palette=cp,legend=True)
	plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
	p.axes.set_title("UMAP plot of cell-mixture scanpy label",fontsize=30)
	p.set_xlabel("UMAP1",fontsize=20)
	p.set_ylabel("UMAP2",fontsize=20)
	plt.savefig(args_sample+'_scanpy_umap.png');plt.close()

	# df_scanpy_pcs = pd.DataFrame(adata.obsm['X_pca'])
	# df_scanpy_pcs['index'] = df_scanpy['index']
	# df_scanpy_pcs = df_scanpy_pcs[['index']+[int(x) for x in df_scanpy_pcs.columns[:-1]]]
	# df_scanpy_pcs.to_pickle(args_sample+'scanpy_raw_pipeline_123k_50PCS_hvariable_genes.pkl')


##################################################################
############ 
##################################################################
