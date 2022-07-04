
def preprocessing():
	import matplotlib.pyplot as plt
	import pandas as pd
	import scanpy as sc

	from pathlib import Path
	server = Path.home().as_posix()
	experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'

	df = pd.read_pickle(server+experiment_home+'data/GSE176078mix/GSE176078mix.counts.FILTERED.pkl')
	experiment_home = server+experiment_home+'output/scanpy/'
	sample='gse_normal'

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
	sc.tl.umap(adata,min_dist=0.1)
	sc.tl.leiden(adata, resolution = 0.2)
	sc.pl.umap(adata, color=['leiden'])
	plt.savefig(experiment_home+ sample+'_scanpy_raw_pipeline_umap.png');plt.close()

	sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
	sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
	plt.savefig(experiment_home+ sample+'_scanpy_raw_pipeline_marker.png');plt.close()

	df_scanpy = adata.to_df()
	df_scanpy = df_scanpy.reset_index()
	
	##save
	df_scanpy.to_csv(experiment_home+ sample+'_scanpy_processed.csv.gz',index=False,compression='gzip')
	
	df_scanpy[['umap1','umap2']] = adata.obsm['X_umap']
	df_scanpy['leiden'] = adata.obs['leiden'].values
	df_scanpy.to_pickle(experiment_home+ sample+'_scanpy_processed.pkl')


def plot_marker():
	import matplotlib.pyplot as plt
	import colorcet as cc
	import seaborn as sns
	plt.rcParams['figure.figsize'] = [15, 10]
	plt.rcParams['figure.autolayout'] = True
	import pandas as pd

	from pathlib import Path
	server = Path.home().as_posix()
	experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'
	experiment_home = server+experiment_home+'output/scanpy/'
	sample='gse_normal'

	
	df = pd.read_pickle(experiment_home+ sample+'_scanpy_processed.pkl')
	df_ann = pd.read_csv(experiment_home+ sample+'_annotated_encode.csv')

	df['labels'] = df_ann['labels']

	# df = df[df['labels'].isin(['B_cell','T_cells','NK_cell','Neuroepithelial_cell','BM'])]
	selected = ['CD8+ T-cells','CD4+ T-cells','B-cells','NK cells','Endothelial cells','Fibroblasts','DC']

	df['labels'] = ['Epithelial cells' if x not in selected else x for x in df['labels']  ]

	df['labels'] = [x.replace('+','').replace(' ','_') for x in df['labels']]
	p = sns.scatterplot(data=df, x='umap1', y='umap2', hue='labels',s=5,legend=True)
	plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
	p.set_xlabel("UMAP1",fontsize=20)
	p.set_ylabel("UMAP2",fontsize=20)
	plt.savefig(experiment_home+ sample+'_scanpy_umap_annotated_encode.png');plt.close()


	fig, ax = plt.subplots(3,4) 
	ax = ax.ravel()


	marker_genes = ['IL7R','CD3D', #CD4
					'CD8A', #CD8
					'LYZ', 'CD14', #CD14
					'MS4A1','CD79A', # B
					'GNLY', 'NKG7',#NK
					'FCER1A',#DC
					'S100A8', 'S100A9'] #myeloid
	
	for i,g in enumerate(marker_genes):
		if g in df.columns:
			print(g)
			h = [ x if x >0 else 0 for x in df[g]]
			sns.scatterplot(data=df, x='umap1', y='umap2', hue=h,s=1,palette="viridis",legend=False,ax=ax[i])
			# ax[i].set_xlabel("UMAP1",fontsize=20)
			# ax[i].set_ylabel("UMAP2",fontsize=20)
			ax[i].set_title(g)
	fig.savefig(experiment_home+ sample+'_scanpy_umap_gene_immune.png');plt.close()

	# b_cells = [ 'CD79A', 'MS4A1']
	# t_cells = ['CD3E', 'CD3G', 'CD3D']
	# t_cell_subsets = ['CD4', 'CD8A', 'CD8B']
	# naive_t_cell = ['SELL', 'CCR7', 'IL7R']
	# myeloid_cells = ['S100A8', 'S100A9', 'CST3']
	# monocytes = ['FCGR3A', 'FCGR3B', 'CD14'] #FCGR3A/B = CD16
	# dendritic_cells = ['FCER1A', 'ITGAX'] #ITGAM = CD11b #ITGAX= CD11c
	# NK_cells = ['NCAM1', 'NKG7', 'CD3G']

	# sc.settings.set_figure_params(dpi=90)
	# #visualize the gene expression as an overlay of the umap
	# #(this way you can visually identify the clusters with a high expression))
	# sc.pl.umap(adata, color = b_cells, color_map = 'viridis', ncols = 3)
	# sc.pl.umap(adata, color = t_cells, color_map = 'viridis', ncols = 3)
	# sc.pl.umap(adata, color = t_cell_subsets, color_map = 'viridis', ncols = 3)
	# sc.pl.umap(adata, color = naive_t_cell, color_map = 'viridis', ncols = 3)
	# sc.pl.umap(adata, color = NK_cells,  color_map = 'viridis',ncols = 3)
	# sc.pl.umap(adata, color = myeloid_cells,  color_map = 'viridis',ncols = 3)
	# sc.pl.umap(adata, color = monocytes, color_map = 'viridis', ncols = 3)
	
	# plt.savefig(experiment_home+ sample+'_scanpy_umap_gene_all.png');plt.close()


	# marker={
	# 	0:'n_T-cells' ,
	# 	1:'n_Epithelial',
	# 	2:'n_Fibroblasts',
	# 	3:'n_Endothelial',
	# 	4:'n_B-cells',
	# 	5:'n_Myo-epithelial',
	# 	6:'n_Epithelial',
	# 	7:'n_Others',
	# 	8:'n_Others',
	# 	9:'n_Others',
	# 	10:'n_Others',
	# 	11:'n_Others'
	# 	}

	# df['leiden'] = df['leiden'].astype(int)
	# h = [ marker[x] for x in df['leiden']]
	# p = sns.scatterplot(data=df, x='umap1', y='umap2', hue=h,s=5,legend=True)
	# plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
	# p.set_xlabel("UMAP1",fontsize=20)
	# p.set_ylabel("UMAP2",fontsize=20)
	# plt.savefig(experiment_home+ sample+'_scanpy_umap_annotated.png');plt.close()

	# df['annot'] = h
	# df[['index','annot']].to_csv(experiment_home+sample+'_annotated.csv.gz',index=False,compression='gzip')


