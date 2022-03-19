import pandas as pd
import matplotlib.pylab as plt
from gen_util.io import read_config
import datatable as dt
from scipy import sparse 
import numpy as np
plt.rcParams["figure.figsize"] = [12.50, 10.50]
plt.rcParams["figure.autolayout"] = True


def select_protein_coding_genes(df):
	params = read_config()
	pcg = params["home"]+params["database"]+params["protein_coding_genes"]
	df_pcg = pd.read_csv(pcg,header=None,names=["gene"])
	drop_columns = [x for x in df.columns[1:-1] if x not in df_pcg["gene"].values]
	df = df.drop(drop_columns,axis=1)
	return df

def scanpy_filter(df):
	
	import scanpy as sc
	
	## convert df to anndata object
	obs = pd.DataFrame(df.iloc[:,0])
	var = pd.DataFrame(index=df.columns[1:-1])
	X = df.iloc[:,1:-1].to_numpy()
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
	plt.savefig("../output/scanpy_etm_tcell_all_pca.png");plt.close()
	
	sc.pl.pca_variance_ratio(adata, log=True)
	plt.savefig("../output/scanpy_tcell_all_pca_var.png");plt.close()

	##optional
	sc.pp.neighbors(adata)
	sc.tl.umap(adata)
	sc.pl.umap(adata, color=['CST3'])
	plt.savefig("../output/scanpy_etm_tcell_all_umap.png");plt.close()

	sc.tl.leiden(adata)
	sc.pl.umap(adata, color=['leiden'])
	plt.savefig("../output/scanpy_tcell_all_umap_leiden.png");plt.close()


	pclust = cellid_to_meta_single(adata.obs["index"].values)

	adata.obs["pclust"]=pclust
	sc.pl.umap(adata, size=2.0, color=['pclust'])
	plt.savefig("../output/scanpy_tcell_all_umap_pclust.png");plt.close()

	marker = ["CD4.c20.Treg.TNFRSF9",\
		"CD4.c06.Tm.ANXA1",\
		"CD4.c01.Tn.TCF7",\
		"CD8.c02.Tm.IL7R",\
		"CD8.c05.Tem.CXCR5",\
		"CD8.c07.Temra.CX3CR1"]

	pclust_mod = [ x if x in marker else "others" for x in pclust]
	adata.obs["pclust_mod"]=pclust_mod
	clrs = ['dodgerblue','purple','green','red','orange','brown','grey']
	sc.pl.umap(adata, size=2.0, color=['pclust_mod'],palette=clrs)
	plt.savefig("../output/scanpy_etm_tcell_all_umap_pclust_selected.png");plt.close()


	pclust_mod2 = [ "CD4" if "CD4." in x else "CD8" for x in pclust]
	adata.obs["pclust_mod2"]=pclust_mod2
	clrs = ['dodgerblue','orange',]
	sc.pl.umap(adata, size=2.0, color=['pclust_mod2'],palette=clrs)
	plt.savefig("../output/scanpy_tcell_all_umap_pclust_selected_CD4_8.png");plt.close()
	# component_check.run_component_analysis(
	#     ,
	#     ,
	#     "tsne",data.split(".")[0])
	# component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"pca",data.split(".")[0])
 
	## get df back from anndata object
	df_scanpy = adata.to_df()
	df_scanpy.columns = adata.var_names
	df_scanpy["cell"] = adata.obs["index"]
	df_scanpy["leiden"] = adata.obs["leiden"]
	df_scanpy["pclust"] = adata.obs["plust"]
	df_scanpy["pclust_mod"] = adata.obs["plust_mod"]

	import umap
	import seaborn as sns
	reducer = umap.UMAP()
	embedding = reducer.fit_transform(df_scanpy.iloc[:,:-2])
	clrs = [sns.color_palette()[int(x)] for x in df_scanpy["leiden"].unique()]
	plt.scatter(
	embedding[:, 0],
	embedding[:, 1],s=0.1)
	plt.gca().set_aspect('equal', 'datalim')
	plt.title('UMAP projection', fontsize=24)
	plt.savefig("../output/scanpy_tcell_all_umap_from_library.png");plt.close()



def filter_minimal(df,mode, cutoff):
	
	print("initial data size : ",df.shape)

	if mode == "tissue":
		# eliminate gene if the total number of counts is < cutoff per tissue type.
		drop_columns =[]
		for tissue in df["sample"].unique():
			drop_columns_sample = [ col for col,val  in df[df["sample"]==tissue].iloc[:,1:-1].sum(axis=0).iteritems() if val < cutoff ]
			for c in drop_columns_sample:
				if c not in drop_columns:
					drop_columns.append(c)
	
	elif mode == "dataset":
		drop_columns = [ col for col,val  in df.iloc[:,1:].sum(axis=0).iteritems() if val < cutoff ]


	df = df.drop(drop_columns,axis=1)
	print("after selecting genes with read counts >= "+str(cutoff)+" cells ",df.shape)
	
	return df


def read_sc_data():
	params = read_config()
	sample_path = params["home"]+params["data"]

	df_combine = pd.DataFrame()
	for sample in params["samples"]:
		print("processing--"+sample)
		dtab = dt.fread(sample_path+sample)
		df = dtab.to_pandas()
		df = df.T 
		df = df.rename(columns=df.iloc[0])
		df = df.iloc[1:].reset_index()
		df = df.rename(columns={"index":"cell"})
		df['sample'] = sample.split('_')[1]+"_"+sample.split('.')[1]
		print(df.head())
		
		## take 500 cells sample per cell typer per tissue type
		df = df.sample(n = 500)
		if sample == params['sample_first']:
			df_combine = df
			print(df.shape,df_combine.shape)
		else:
			df_combine = pd.concat([df_combine, df], axis=0, ignore_index=True)
			print(df.shape,df_combine.shape)
	
	# return df_combine
	# df_combine = df_combine.fillna(0) ###super slow
	df_combine.values[df_combine.isna()] = 0
	df_combine.to_csv("../output/cd4_cd8_500cells_per_tissue_counts.txt.gz",index=False,sep="\t",compression="gzip")

def read_sc_data_sparse(maxcellcount_pertissue=1e5,mingene_totalcount=10):

	params = read_config()
	sample_path = params["home"]+params["data"]

	df_combine = pd.DataFrame()
	for sample in params["samples"]:
		
		print("processing--"+sample)
		
		dtab = dt.fread(sample_path+sample,verbose=False)
		df = dtab.to_pandas() 
		df = df.T
		df = df.rename(columns=df.iloc[0])
		df = df.iloc[1:].reset_index()
		
		df = filter_minimal(df,"dataset", mingene_totalcount)

		if df.shape[0]>maxcellcount_pertissue:
			df = df.sample(n=maxcellcount_pertissue)
		if sample == params['sample_first']:
			df_combine = df
			print(df.shape,df_combine.shape)
		else:
			print("concatinating...")
			df_combine = pd.concat([df_combine, df], axis=0, ignore_index=True)
			print(df.shape,df_combine.shape)
	del df
	
	print("fill nans as zero...")
	df_combine.values[df_combine.isna()] = 0
	
	print("processing--creating coo sparse matrix file")

	label = df_combine.iloc[:,0].to_numpy()
	np.savez_compressed(params["home"]+params["data"]+params["sparse_label"]+"_"+str(maxcellcount_pertissue), idx=label,allow_pickle=True)

	df_combine = df_combine.iloc[:,1:]
	immune_signatures = get_immune_genes()
	immune_index = [(x,y) for x,y in enumerate(df_combine.columns) if y in immune_signatures]
	np.savez_compressed(params["home"]+params["data"]+params["sparse_label"]+"_immune_index_"+str(maxcellcount_pertissue), idx=immune_index, allow_pickle=True)

	S = sparse.coo_matrix(df_combine.to_numpy())
	idx, idy, val = sparse.find(S)
	d = df_combine.shape
	np.savez_compressed(params["home"]+params["data"]+params["sparse_data"]+"_"+str(maxcellcount_pertissue), idx=idx,idy=idy,val=val,shape=d,allow_pickle=True)

def select_cells():
	'''
	Create a test data for algorithm development from Zheng et al paper (GSE156728) 
	selected data - tissue type, cell type, cluster.
	'''
	
	celltype_cluster = ['CD4_c18','CD4_c13']
	tissue_type = ['BC', 'BCL', 'ESCA', 'MM', 'PACA', 'RC', 'THCA', 'UCEC']
	
	params = read_config()
	meta_path = params["home"]+params["database"]
	
	df = pd.read_csv(meta_path+params["metadata"],sep="\t")
	
	df["meta.cluster.sf"] =[x.split('.')[0]+"_"+x.split('.')[1] for x in df["meta.cluster"]]
	df = df[(df["meta.cluster.sf"].isin(celltype_cluster))]
	df = df[df["cancerType"].isin(tissue_type)]
	df['cluster'] = [x+"_"+y for x,y in zip(df["cancerType"],df["meta.cluster.sf"])]
	df = df[["cellID","cluster"]]
	df.to_csv("../output/all_tissue_cd4_c13_c18_cell_barcodes.csv",sep="\t",index=False)

def create_test_data():
	## prepare data matrix for development
	df = read_sc_data()
	df = df.fillna(0)
	## read selected cell barcodes
	df_bc = pd.read_csv("../output/all_tissue_cd4_c13_c18_cell_barcodes.csv",sep="\t")
	df_selected = pd.merge(df_bc,df,on="cell",how="left") 
	'''IMPORTANT:memory reaches 90%, 80k by 30k merge with 5k by 2'''
	df_selected= df_selected[ [ col for col in df_selected.columns if col != "cluster" ] + ["cluster"]]
	df_selected.to_csv("../output/cd4_c13_c18_counts.txt.gz",index=False,sep="\t",compression="gzip")


def cellid_to_meta(cell_pairs,args):
	meta_path = args.home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep="\t")

	cell_pairs_meta = []
	for cp in cell_pairs:
		types = df_meta.loc[df_meta["cellID"].isin(cp),["meta.cluster"]].values
		cell_pairs_meta.append(types.flatten()[0]+"/"+types.flatten()[1])

	return cell_pairs_meta

def cellid_to_meta_single(cell_ids,args):
	meta_path = args.home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep="\t")
	df_meta = df_meta[["cellID","meta.cluster"]]

	df = pd.DataFrame(cell_ids,columns=["cellID"])

	dfjoin = pd.merge(df,df_meta,on="cellID",how="left")

	return dfjoin["meta.cluster"].values

def get_immune_genes():
	params = read_config()
	meta_path = params["home"]+params["database"]
	df_meta = pd.read_csv(meta_path+params["immune_signature_genes"],sep="\t")
	return df_meta["signature_genes"].values

def create_lr_exp_data(args):
	lr_db=args.home+args.database+args.lr_db
	df = pd.read_csv(lr_db,sep="\t")
	receptors = list(df.receptor_gene_symbol.unique())
	ligands = list(df.ligand_gene_symbol.unique())

	df_exp=pd.read_pickle(args.home+args.input+args.raw_data)

	fname = args.home+args.input+args.raw_data
	l_fname = fname.replace(".pkl",".ligands.pkl")
	r_fname = fname.replace(".pkl",".receptors.pkl")

	df_exp[["index"]+[ x for x in ligands if x in df_exp.columns ]].to_pickle(l_fname)
	df_exp[["index"]+[ x for x in receptors if x in df_exp.columns ]].to_pickle(r_fname)

def create_lr_mat(args):
	lr_db=args.home+args.database+args.lr_db
	df = pd.read_csv(lr_db,sep="\t")
	dflrmat = df.groupby(['ligand_gene_symbol','receptor_gene_symbol']).agg(['count'])['lr_pair']
	dflrmat = dflrmat.unstack(fill_value=0)
	dflrmat.columns = dflrmat.columns.droplevel(0)
	fname = args.home+args.input+args.raw_data
	lr_fname = fname.replace(".pkl",".ligands_receptors_mat_815_780.pkl")
	dflrmat.to_pickle(lr_fname)


