import os
import pandas as pd
import datatable as dt
from scipy import sparse 
import numpy as np
import logging
logger = logging.getLogger(__name__)


def select_protein_coding_genes(args,df):
	df_pcg = pd.read_csv(os.environ['args_home']+args.database+args.protein_coding_genes,header=None,names=['gene'])
	drop_columns = [x for x in df.columns[1:-1] if x not in df_pcg['gene'].values]
	df = df.drop(drop_columns,axis=1)
	return df

def get_ligand_receptor_genes(args):
	df = pd.read_csv( os.environ['args_home']+args.database+args.lr_db,sep='\t')
	return list(df.receptor_gene_symbol.unique())+list(df.ligand_gene_symbol.unique())

def filter_minimal(df,cutoff):
	
	logger.info('initial data size : '+str(df.shape[0])+'x'+str(df.shape[1]))

	drop_columns = [ col for col,val  in df.iloc[:,1:].sum(axis=0).iteritems() if val < cutoff ]

	logger.info('genes to filter based on mincout cutoff - '+ str(len(drop_columns)))

	for mt_g in [x for x in df.columns if 'MT-' in x]:
		drop_columns.append(mt_g)

	logger.info('adding mitochondrial genes - '+ str(len(drop_columns)))

	for spk_g in [x for x in df.columns if 'ERCC' in x]:
		drop_columns.append(spk_g)

	logger.info('adding spikes - '+ str(len(drop_columns)))

	df = df.drop(drop_columns,axis=1)

	logger.info('after all minimal filtering..'+str(df.shape[0])+'x'+str(df.shape[1]))

	return df

def get_immune_genes(args):
	args_home = os.environ['args_home']
	df_meta = pd.read_csv(args_home+args.database+args.immune_signature_genes,sep="\t")

	df_meta = df_meta[
				( df_meta['species'].isin(['Mm Hs','Hs']) ) &
			    ( df_meta['organ'].isin(['Immune system']) )
			]
			
	return df_meta['official gene symbol'].unique()

def tcell_signature_genes(args_home,args):
	df_meta = pd.read_csv(args_home+args.database+'pancancer_tcell_signature_gene_all_clusters.tsv')
	return df_meta['signature_genes'].values

def lr_exp_data(args,df_exp):
	args_home = os.environ['args_home']
	df = pd.read_csv( args_home + args.database+args.lr_db,sep='\t')
	receptors = list(df.receptor_gene_symbol.unique())
	ligands = list(df.ligand_gene_symbol.unique())

	df_exp[['index']+[ x for x in ligands if x in df_exp.columns ]].to_pickle(
		args_home+args.input+args.raw_l_data)
	df_exp[['index']+[ x for x in receptors if x in df_exp.columns ]].to_pickle(
		args_home+args.input+args.raw_r_data)

	pd.DataFrame([ x for x in ligands if x in df_exp.columns ]).to_pickle(args_home+args.input+args.raw_l_data_genes)
	pd.DataFrame([ x for x in receptors if x in df_exp.columns ]).to_pickle(args_home+args.input+args.raw_r_data_genes)

def create_lr_mat(args):
	df = pd.read_csv( os.environ['args_home']+args.database+args.lr_db,sep='\t')
	dflrmat = df.groupby(['ligand_gene_symbol','receptor_gene_symbol']).agg(['count'])['lr_pair']
	dflrmat = dflrmat.unstack(fill_value=0)
	dflrmat.columns = dflrmat.columns.droplevel(0)
	fname = os.environ['args_home']+args.input+args.raw_lr_data
	dflrmat.to_pickle(fname)

