import os
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np

def preprocessing():

	from pathlib import Path
	server = Path.home().as_posix()
	experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'
	experiment_home = server+experiment_home
	sample='gse_normal'

	df = pd.read_pickle(experiment_home+'data/GSE176078mix/GSE176078mix.counts.FILTERED.pkl')
	df['label'] = [x.split('_')[len(x.split('_'))-1] for x in df['index']]
	df = df[df['label']=='GSE164898']
	df = df.drop(columns=['label'])

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
	plt.savefig(experiment_home+ sample+'_scanpy_raw_pipeline_pca.png');plt.close()
	
	sc.pl.pca_variance_ratio(adata, n_pcs=50,log=True)
	plt.savefig(experiment_home+ sample+'_scanpy_raw_pipeline_pca_var.png');plt.close()

	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
	sc.tl.umap(adata)
	sc.tl.leiden(adata, resolution = 0.1)
	sc.pl.umap(adata, color=['leiden'])
	plt.savefig(experiment_home+ sample+'_scanpy_raw_pipeline_umap.png');plt.close()

	sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
	sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
	plt.savefig(experiment_home+ sample+'_scanpy_raw_pipeline_marker.png');plt.close()

	df_scanpy = adata.to_df()
	df_scanpy = df_scanpy.reset_index()
	df_scanpy[['umap1','umap2']] = adata.obsm['X_umap']
	df_scanpy['leiden'] = adata.obs['leiden'].values
	df_scanpy.to_pickle(experiment_home+ sample+'_scanpy_processed.pkl')


def plot_marker():
	
	import colorcet as cc
	import seaborn as sns
	plt.rcParams['figure.figsize'] = [15, 10]
	plt.rcParams['figure.autolayout'] = True


	from pathlib import Path
	server = Path.home().as_posix()
	experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'
	experiment_home = server+experiment_home
	sample='gse_normal'

	
	df = pd.read_pickle(experiment_home+ sample+'_scanpy_processed.pkl')

	cp = sns.color_palette(cc.glasbey, n_colors=len(df['leiden'].unique()))
	p = sns.scatterplot(data=df, x='umap1', y='umap2', hue='leiden',s=1,palette=cp,legend=False)
	p.set_xlabel("UMAP1",fontsize=20)
	p.set_ylabel("UMAP2",fontsize=20)
	plt.savefig(experiment_home+ sample+'_scanpy_umap_t.png');plt.close()


	h = [ x if x >0 else 0 for x in df['NR4A2']]
	p = sns.scatterplot(data=df, x='umap1', y='umap2', hue=h,s=1,legend=False)
	p.set_xlabel("UMAP1",fontsize=20)
	p.set_ylabel("UMAP2",fontsize=20)
	plt.savefig(experiment_home+ sample+'_scanpy_umap_gene.png');plt.close()


	marker={
		0:'n_T-cells' ,
		1:'n_Epithelial',
		2:'n_Fibroblasts',
		3:'n_Endothelial',
		4:'n_B-cells',
		5:'n_Myo-epithelial',
		6:'n_Epithelial',
		7:'n_Others',
		8:'n_Others',
		9:'n_Others',
		10:'n_Others',
		11:'n_Others'
		}

	df['leiden'] = df['leiden'].astype(int)
	h = [ marker[x] for x in df['leiden']]
	p = sns.scatterplot(data=df, x='umap1', y='umap2', hue=h,s=5,legend=True)
	plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
	p.set_xlabel("UMAP1",fontsize=20)
	p.set_ylabel("UMAP2",fontsize=20)
	plt.savefig(experiment_home+ sample+'_scanpy_umap_annotated.png');plt.close()

	df['annot'] = h
	df[['index','annot']].to_csv(experiment_home+sample+'_annotated.csv.gz',index=False,compression='gzip')


