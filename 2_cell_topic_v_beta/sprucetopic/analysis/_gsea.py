import pandas as pd
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