import os
import pandas as pd
import numpy as np
from premodel import prep_tcells as prep

def topic_top_genes(args,top_n=5):

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

	top_genes = []
	top_genes = prep.generate_gene_vals(df_beta1,top_n,top_genes,'T cell')
	top_genes = prep.generate_gene_vals(df_beta2,top_n,top_genes,'non T cell')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])

	db = args_home+args.database+args.tcell_signature_genes

	df_cd4 = pd.read_excel(db,sheet_name='CD4',header=1,usecols='A,B,C' )
	df_cd4 = df_cd4.sort_values('comb.ES',ascending=False).drop_duplicates('geneSymbol')

	df_cd8 = pd.read_excel(db,sheet_name='CD8',header=1,usecols='A,B,C' )
	df_cd8 = df_cd8.sort_values('comb.ES',ascending=False).drop_duplicates('geneSymbol')

	df_cd48_top_genes = pd.DataFrame(df_top_genes['Gene'].unique(),columns=['Gene'])

	df_cd48_top_genes = pd.merge(df_cd48_top_genes,df_cd4[['geneSymbol','cluster.name']],how='left',right_on='geneSymbol',left_on='Gene')
	df_cd48_top_genes = pd.merge(df_cd48_top_genes,df_cd8[['geneSymbol','cluster.name']],how='left',right_on='geneSymbol',left_on='Gene')

	df_cd48_top_genes = df_cd48_top_genes[['Gene', 'cluster.name_x', 'cluster.name_y']]

	df_cd48_top_genes.columns=['Gene','CD4','CD8']
	df_cd48_top_genes['CD4'] = df_cd48_top_genes['CD4'].fillna('NotMarker')
	df_cd48_top_genes['CD8'] = df_cd48_top_genes['CD8'].fillna('NotMarker')

	df_cd48_top_genes.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'top_'+str(top_n)+'_genes_topic_CDMarkers.tsv',sep='\t',index=False)

	df_top_genes.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'top_'+str(top_n)+'_genes_topic.tsv',sep='\t',index=False)


def tcell_sample_cells_with_latent(args,cell_n=100):

	args_home = os.environ['args_home']

	df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_hh_data.tsv',sep='\t',compression='gzip')


	df_h['cluster'] = prep.cellid_to_meta_single(df_h['cell'].values,args)

	markers=['CD4.c04.Tn.il7r', \
				'CD4.c07.Tm.ANXA2', \
				'CD4.c12.Tem.GZMK', \
				'CD4.c13.Temra.CX3CR1', \
				'CD4.c14.Th17.SLC4A10', \
				'CD4.c16.Tfh.CXCR5', \
				'CD4.c17.TfhTh1.CXCL13',\
				'CD4.c20.Treg.TNFRSF9',\
				'CD4.c21.Treg.OAS1',\
				'CD4.c22.ISG.IFIT1',\
				'CD8.c01.Tn.MAL',\
				'CD8.c05.Tem.CXCR5',\
				'CD8.c06.Tem.GZMK',\
				'CD8.c07.Temra.CX3CR1',\
				'CD8.c08.Tk.TYROBP',\
				'CD8.c11.Tex.PDCD1',\
				'CD8.c15.ISG.IFIT1',\
				'CD8.c16.MAIT.SLC4A10',\
				'CD8.c17.Tm.NME1'
			   ]

	df_h = df_h[df_h['cluster'].isin(markers)]

	df_h_sample = df_h.groupby('cluster').sample(cell_n, random_state=1)

	df_h_sample.columns = [ x.replace('hh','k') for x in df_h_sample.columns]
	df_h_sample.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'hh_cell_topic_sample_CD4CD8.tsv',sep='\t',index=False)

def top_genes_marker_correlation(args,top_n=5):

	args_home = os.environ['args_home']

	df_genes=pd.read_pickle(args_home+args.input+args.raw_data_genes)

	immuneindex_file = args_home+args.input+args.nbr_model['sparse_immune_label']
	immunedat = np.load(immuneindex_file,allow_pickle=True)
	immuneindx = np.array([x[0] for x in immunedat['idx']]).astype(np.int32)
	otherindx = np.array([x for x in range(df_genes.shape[0]) if x not in immuneindx]).astype(np.int32)

	df_beta1 = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_beta1_data.tsv',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_beta2_data.tsv',sep='\t',compression='gzip')

	df_beta1.columns = df_genes.iloc[immuneindx][0].values
	df_beta2.columns = df_genes.iloc[otherindx][0].values

	db = args_home+args.database+args.tcell_signature_genes

	tag='CD8'
	df_cd4 = pd.read_excel(db,sheet_name=tag,header=1,usecols='A,B,C' )
	df_cd4 = df_cd4.sort_values('comb.ES',ascending=False)

	df_cd4_top5 = df_cd4.groupby('cluster.name').head(5).reset_index(drop=True)
	df_cd4_top5 = df_cd4_top5.sort_values('cluster.name')

	topic_clust_query = []
	clusters = df_cd4_top5['cluster.name'].unique()
	topic_clust = []
	for topic in range(df_beta1.shape[0]):
		topic_genes = df_beta1.T.iloc[:,topic].sort_values(0,ascending=False)[:top_n].index

		topic_genes_vals = df_beta1.T.iloc[:,topic].sort_values(0,ascending=False)[:top_n].values
		cell_clust = []
		for clust in clusters :
			clust_genes = df_cd4_top5[df_cd4_top5['cluster.name']==clust]['geneSymbol'].values
			clust_genes_vals = df_beta1.loc[topic,clust_genes].values
			cor = np.corrcoef(topic_genes_vals,clust_genes_vals)[0,1]
			cell_clust.append(cor)
			topic_clust_query.append(
				[topic,
				topic_genes,
				topic_genes_vals,
				clust,
				clust_genes,
				clust_genes_vals,
				cor]
			)
		topic_clust.append(cell_clust)

	df_topic_cluster_corr = pd.DataFrame(topic_clust)
	df_topic_cluster_corr.columns = clusters

	df_topic_cluster_corr.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'top_'+str(top_n)+'_topic_cluster_corr'+tag+'.tsv',sep='\t',index=False)

	pd.DataFrame(topic_clust_query).to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'top_'+str(top_n)+'_topic_cluster_corr_query'+tag+'.tsv',sep='\t',index=False)

def get_lr_pair_topic_cd48_top_genes(args,top_n=1):

	args_home = os.environ['args_home']
	l_fname = args_home+args.input+args.raw_l_data_genes
	r_fname = args_home+args.input+args.raw_r_data_genes

	ligands = pd.read_pickle(l_fname)[0].values
	receptors = pd.read_pickle(r_fname)[0].values

	df_beta1 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta1_data.tsv',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta2_data.tsv',sep='\t',compression='gzip')

	df_beta1.columns = receptors
	df_beta2.columns = ligands

	top_genes = []
	top_genes = prep.generate_gene_vals(df_beta1,top_n,top_genes,'receptors')
	top_genes = prep.generate_gene_vals(df_beta2,top_n,top_genes,'ligands')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])
	top_genes = df_top_genes['Gene'].unique()

	db = args_home+args.database+args.tcell_signature_genes

	df_cd4 = pd.read_excel(db,sheet_name='CD4',header=1,usecols='A,B,C' )
	df_cd4 = df_cd4.sort_values('comb.ES',ascending=False).drop_duplicates('geneSymbol')

	df_cd8 = pd.read_excel(db,sheet_name='CD8',header=1,usecols='A,B,C' )
	df_cd8 = df_cd8.sort_values('comb.ES',ascending=False).drop_duplicates('geneSymbol')

	cd48_top_genes = list(df_cd4.geneSymbol.unique())+list(df_cd8.geneSymbol.unique())

	topic_cd48_top_genes = [ x for x in cd48_top_genes if x in top_genes]

	df_db = pd.read_csv( args_home + args.database+args.lr_db,sep='\t', usecols=['lr_pair'])

	lr_topic = []
	lr_pair = []
	for lr_p in df_db['lr_pair']:
		l,r = lr_p.split('_')[0],lr_p.split('_')[1]
		if l in df_beta2.columns and r in df_beta1.columns and \
			l in topic_cd48_top_genes and r in topic_cd48_top_genes:
			lr_topic.append((df_beta2[l]+df_beta1[r])/2)
			lr_pair.append(lr_p)
	df_lr_topic = pd.DataFrame(lr_topic)
	df_lr_topic.index=lr_pair
	df_lr_topic.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'top_'+str(top_n)+'_lrpair_topic_cd4cd8_top_genes.tsv',sep='\t')

