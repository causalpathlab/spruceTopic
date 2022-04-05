import os
import pandas as pd
import matplotlib.pylab as plt
from gen_util.io import read_config
import datatable as dt
from scipy import sparse 
import numpy as np
plt.rcParams['figure.figsize'] = [12.50, 10.50]
plt.rcParams['figure.autolayout'] = True


def select_protein_coding_genes(df):
	params = read_config()
	pcg = params['home']+params['database']+params['protein_coding_genes']
	df_pcg = pd.read_csv(pcg,header=None,names=['gene'])
	drop_columns = [x for x in df.columns[1:-1] if x not in df_pcg['gene'].values]
	df = df.drop(drop_columns,axis=1)
	return df


def get_ligand_receptor_genes(args):
	df = pd.read_csv( os.environ['args_home']+args.database+args.lr_db,sep='\t')
	return list(df.receptor_gene_symbol.unique())+list(df.ligand_gene_symbol.unique())

def filter_minimal(args,df,mode, cutoff):
	
	print('initial data size : ',df.shape)

	ligand_receptors = get_ligand_receptor_genes(args)
	drop_columns =[]

	if mode == 'tissue':
		# eliminate gene if the total number of counts is < cutoff per tissue type.
		for tissue in df['sample'].unique():
			drop_columns_sample = [ col for col,val  in df[df['sample']==tissue].iloc[:,1:-1].sum(axis=0).iteritems() if val < cutoff ]
			for c in drop_columns_sample:
				if c not in drop_columns:
					drop_columns.append(c)
	
	elif mode == 'dataset':
		drop_columns = [ col for col,val  in df.iloc[:,1:].sum(axis=0).iteritems() if val < cutoff ]

	print('genes to filter based on mincout cutoff - ', len(drop_columns))

	for mt_g in [x for x in df.columns if 'MT-' in x]:
		drop_columns.append(mt_g)

	print('adding mitochondrial genes - ' ,len(drop_columns))

	for spk_g in [x for x in df.columns if 'ERCC' in x]:
		drop_columns.append(mt_g)

	print('adding spikes - '  ,len(drop_columns))

	for lr_g in [x for x in df.columns if x in ligand_receptors]:
		drop_columns.append(lr_g)

	print('adding ligand-receptors - '  ,len(drop_columns))


	df = df.drop(drop_columns,axis=1)
	print('after all minimal filtering..',df.shape)


	return df

def read_sc_data():
	params = read_config()
	sample_path = params['home']+params['data']

	df_combine = pd.DataFrame()
	for sample in params['samples']:
		print('processing--'+sample)
		dtab = dt.fread(sample_path+sample)
		df = dtab.to_pandas()
		df = df.T 
		df = df.rename(columns=df.iloc[0])
		df = df.iloc[1:].reset_index()
		df = df.rename(columns={'index':'cell'})
		df['sample'] = sample.split('_')[1]+'_'+sample.split('.')[1]
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
	df_combine.to_csv('../output/cd4_cd8_500cells_per_tissue_counts.txt.gz',index=False,sep='\t',compression='gzip')

def read_sc_data_sparse(args,maxcellcount_pertissue=1e5,mingene_totalcount=10):
	args_home = os.environ['args_home'] 
	df_combine = pd.DataFrame()
	sample_path = args_home+args.input+args.nbr_model['in']
	sample_procpath = args_home+args.input

	for sample in args.samples[0:3]:
		
		print('processing--'+sample)
		
		dtab = dt.fread(sample_path+sample,verbose=False)
		df = dtab.to_pandas() 
		df = df.T
		df = df.rename(columns=df.iloc[0])
		df = df.iloc[1:].reset_index()
		
		df = filter_minimal(args,df,'dataset', mingene_totalcount)

		df = df.sample(n=10)
		
		if df.shape[0]>maxcellcount_pertissue:
			df = df.sample(n=maxcellcount_pertissue)
		if sample == args.sample_first:
			df_combine = df
			print(df.shape,df_combine.shape)
		else:
			print('concatinating...')
			df_combine = pd.concat([df_combine, df], axis=0, ignore_index=True)
			print(df.shape,df_combine.shape)
	del df
	print('fill nans as zero...')
	df_combine.values[df_combine.isna()] = 0

	## save data 
	df_combine.to_pickle(sample_procpath+args.raw_data)
	pd.DataFrame(df_combine.columns[1:]).to_pickle(sample_procpath+args.raw_data_genes)

	print('processing--creating coo sparse matrix file')

	label = df_combine.iloc[:,0].to_numpy()
	np.savez_compressed(sample_procpath+args.nbr_model['sparse_label'], idx=label,allow_pickle=True)

	df_combine = df_combine.iloc[:,1:]
	immune_signatures = get_tcell_signature_genes(args_home,args)
	immune_index = [(x,y) for x,y in enumerate(df_combine.columns) if y in immune_signatures]
	np.savez_compressed(sample_procpath + args.nbr_model['sparse_immune_label'], idx=immune_index, allow_pickle=True)

	S = sparse.coo_matrix(df_combine.to_numpy())
	idx, idy, val = sparse.find(S)
	d = df_combine.shape
	np.savez_compressed(sample_procpath + args.nbr_model['sparse_data'], idx=idx,idy=idy,val=val,shape=d,allow_pickle=True)

def select_cells():
	'''
	Create a test data for algorithm development from Zheng et al paper (GSE156728) 
	selected data - tissue type, cell type, cluster.
	'''
	
	celltype_cluster = ['CD4_c18','CD4_c13']
	tissue_type = ['BC', 'BCL', 'ESCA', 'MM', 'PACA', 'RC', 'THCA', 'UCEC']
	
	params = read_config()
	meta_path = params['home']+params['database']
	
	df = pd.read_csv(meta_path+params['metadata'],sep='\t')
	
	df['meta.cluster.sf'] =[x.split('.')[0]+'_'+x.split('.')[1] for x in df['meta.cluster']]
	df = df[(df['meta.cluster.sf'].isin(celltype_cluster))]
	df = df[df['cancerType'].isin(tissue_type)]
	df['cluster'] = [x+'_'+y for x,y in zip(df['cancerType'],df['meta.cluster.sf'])]
	df = df[['cellID','cluster']]
	df.to_csv('../output/all_tissue_cd4_c13_c18_cell_barcodes.csv',sep='\t',index=False)

def cellid_to_meta(cell_pairs,args):
	meta_path = args.home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep='\t')

	cell_pairs_meta = []
	for cp in cell_pairs:
		types = df_meta.loc[df_meta['cellID'].isin(cp),['meta.cluster']].values
		cell_pairs_meta.append(types.flatten()[0]+'/'+types.flatten()[1])

	return cell_pairs_meta

def cellid_to_meta_single(cell_ids,args):
	args_home = os.environ['args_home'] 
	meta_path = args_home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep='\t')
	df_meta = df_meta[['cellID','meta.cluster']]

	df = pd.DataFrame(cell_ids,columns=['cellID'])

	dfjoin = pd.merge(df,df_meta,on='cellID',how='left')

	return dfjoin['meta.cluster'].values

def get_immune_genes():
	params = read_config()
	meta_path = params['home']+params['database']
	df_meta = pd.read_csv(meta_path+params['immune_signature_genes'],sep='\t')
	return df_meta['signature_genes'].values

def get_tcell_signature_genes(args_home,args):
	df_meta = pd.read_csv(args_home+args.database+'pancancer_tcell_signature_gene_all_clusters.tsv')
	return df_meta['signature_genes'].values

def create_lr_exp_data(args):
	lr_db=args.home+args.database+args.lr_db
	df = pd.read_csv(lr_db,sep='\t')
	receptors = list(df.receptor_gene_symbol.unique())
	ligands = list(df.ligand_gene_symbol.unique())

	df_exp=pd.read_pickle(args.home+args.input+args.raw_data)

	fname = args.home+args.input+args.raw_data
	l_fname = fname.replace('.pkl','.ligands.pkl')
	r_fname = fname.replace('.pkl','.receptors.pkl')

	df_exp[['index']+[ x for x in ligands if x in df_exp.columns ]].to_pickle(l_fname)
	df_exp[['index']+[ x for x in receptors if x in df_exp.columns ]].to_pickle(r_fname)

def create_lr_mat(args):
	lr_db=args.home+args.database+args.lr_db
	df = pd.read_csv(lr_db,sep='\t')
	dflrmat = df.groupby(['ligand_gene_symbol','receptor_gene_symbol']).agg(['count'])['lr_pair']
	dflrmat = dflrmat.unstack(fill_value=0)
	dflrmat.columns = dflrmat.columns.droplevel(0)
	fname = args.home+args.input+args.raw_data
	lr_fname = fname.replace('.pkl','.ligands_receptors_mat_815_780.pkl')
	dflrmat.to_pickle(lr_fname)

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