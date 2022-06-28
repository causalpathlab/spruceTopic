import pandas as pd
import gseapy as gp
import numpy as np
from scipy.stats import hypergeom


def beta_z_score(df_mean,df_sd,smean=0.0,n=1e-3):
	return (df_mean-smean)/(df_sd/np.sqrt(n))

def get_degs(genes,dfz,cis = 2.575):
	return [genes[ [ y for y,x in enumerate(dfz.iloc[t,:].values) if abs(x) > cis ] ] for t in dfz.index ]

def hypergeom_test(sp,df_gset):

	dfz = beta_z_score(sp.cell_topic.beta_mean,np.sqrt(np.exp(sp.cell_topic.beta_var)))
	topic_degs = get_degs(sp.data.raw_data_genes,dfz)

	N=len(sp.data.raw_data_genes)
	gse = []
	for pathway in df_gset['gs_name'].unique():
		pathway_genes = df_gset[df_gset['gs_name']==pathway]['gene_symbol'].values
		m = len(pathway_genes)
		for ti,t_deg in enumerate(topic_degs):
			k = len(t_deg)
			q = len([x for x in t_deg if x in pathway_genes])
			pval = 1- hypergeom.cdf(q,N,m,k)
			gse.append([pathway,ti,m,k,q,pval])
	return pd.DataFrame(gse,columns=['pathway','topic','m','k','q','pval'])

def gsea(spr):
	# select gene set database
	# https://maayanlab.cloud/Enrichr/#libraries
	df = spr.cell_topic.beta_mean.copy()

	for i in range(df.shape[0]):
		rnk = pd.DataFrame(df.iloc[i,:].sort_values(axis=0,ascending=False)).reset_index().rename(columns={'index':0,0:1}).iloc[0:500,]
		res = gp.prerank(rnk=rnk,
						gene_sets='CellMarker_Augmented_2021',
						processes=4,
						permutation_num=100,
						outdir=None,
						no_plot=True,seed=6).res2d.sort_values('fdr').iloc[0:2,3].reset_index()
		print(i)
		print(res)