import os
import pandas as pd
import numpy as np

def generate_gene_vals(df,top_n,top_genes,label):

	top_genes_collection = []
	for x in range(df.shape[0]):
		gtab = df.T.iloc[:,x].sort_values(0,ascending=False)[:top_n].reset_index()
		gtab.columns = ['gene','val']
		genes = gtab['gene'].values
		for g in genes:
			if g not in top_genes_collection:
				top_genes_collection.append(g)

	for g in top_genes_collection:
		for i,x in enumerate(df[g].values):
			top_genes.append(['k'+str(i),label,'g'+str(i+1),g,x])

	return top_genes

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
	top_genes = generate_gene_vals(df_beta1,top_n,top_genes,'immune')
	top_genes = generate_gene_vals(df_beta2,top_n,top_genes,'non-immune')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])

	df_top_genes.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'top_'+str(top_n)+'_genes_topic.tsv',sep='\t',index=False)


def get_top_lr_genes(args,top_n=5):

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
	top_genes = generate_gene_vals(df_beta1,top_n,top_genes,'receptors')
	top_genes = generate_gene_vals(df_beta2,top_n,top_genes,'ligands')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])
	df_top_genes.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'top_'+str(top_n)+'_genes_topic.tsv',sep='\t',index=False)

def get_lr_pair_topic(args,top_n=5):

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
	top_genes = generate_gene_vals(df_beta1,top_n,top_genes,'receptors')
	top_genes = generate_gene_vals(df_beta2,top_n,top_genes,'ligands')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])
	top_genes = df_top_genes['Gene'].unique()

	df_db = pd.read_csv( args_home + args.database+args.lr_db,sep='\t', usecols=['lr_pair'])

	lr_topic = []
	lr_pair = []
	for lr_p in df_db['lr_pair']:
		l,r = lr_p.split('_')[0],lr_p.split('_')[1]
		if l in df_beta2.columns and r in df_beta1.columns and \
			l in top_genes and r in top_genes:
			lr_topic.append((df_beta2[l]+df_beta1[r])/2)
			lr_pair.append(lr_p)
	df_lr_topic = pd.DataFrame(lr_topic)
	df_lr_topic.index=lr_pair
	df_lr_topic.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'top_'+str(top_n)+'_lrpair_topic.tsv',sep='\t')

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
	top_genes = generate_gene_vals(df_beta1,top_n,top_genes,'receptors')
	top_genes = generate_gene_vals(df_beta2,top_n,top_genes,'ligands')

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





