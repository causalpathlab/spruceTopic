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

def gse_celltopic(spr):
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


def gse_interactiontopic(spr,df_db):
	# select gene set database
	# https://maayanlab.cloud/Enrichr/#libraries
	df_r = spr.interaction_topic.beta_l
	df_l = spr.interaction_topic.beta_r


	lr_topic = []
	lr_pair = []
	for lr_p in df_db['lr_pair']:
		l,r = lr_p.split('_')[0],lr_p.split('_')[1]
		if l in df_l.columns and r in df_r.columns:
			lr_topic.append((df_l.loc[:,l].values + df_r.loc[:,r].values)/2)
			lr_pair.append(lr_p)
	df_lr_topic = pd.DataFrame(lr_topic)
	df_lr_topic.index=lr_pair
	
	df = df_lr_topic.T

	df_grps = pd.DataFrame()

	for i in range(df.shape[0]):
		rnk = pd.DataFrame(df.iloc[i,:].sort_values(axis=0,ascending=False)).reset_index().rename(columns={'index':0,0:1})

		rnk.columns = [0,1]
		# rnk = rnk[rnk[1]>0.0]

		print(rnk.shape)


		rnk['r'] = [x.split('_')[0] for x in rnk[0]]
		rnk['l'] = [x.split('_')[1] for x in rnk[0]]

		
		rnk_combined = pd.concat( [ rnk[['r',1]].rename(columns={'r':'gene'}) , rnk[['l',1]].rename(columns={'l':'gene'}) ],axis=0, ignore_index=True)

		rnk_combined = rnk_combined.sort_values(1,ascending=False).rename(columns={'gene':0}).reset_index(drop=True)
		
		res = gp.prerank(rnk=rnk_combined,
						gene_sets='MSigDB_Hallmark_2020',
						processes=4,
						permutation_num=100,
						outdir=None,
						no_plot=True,seed=6).res2d.sort_values('fdr').iloc[0:5,3].reset_index()
		
		res = res[res['fdr']<0.01]
		
		if res.shape[0]>0:
			res['topic']=i
			df_grps = pd.concat([df_grps, res], axis=0, ignore_index=True)

	return df_grps

def gse_interactiontopic_v2(spr,df_db):
	# select gene set database
	# https://maayanlab.cloud/Enrichr/#libraries
	df_r = spr.interaction_topic.beta_l
	df_l = spr.interaction_topic.beta_r

	
	df = pd.concat([df_r,df_l],axis=1)

	df_grps = pd.DataFrame()

	for i in range(df.shape[0]):
		# rnk = pd.DataFrame(df.iloc[i,:].sort_values(axis=0,ascending=False)).reset_index().rename(columns={'index':0,0:1})
		rnk_combined = pd.DataFrame(df.iloc[i,:].sort_values(axis=0,ascending=False)).reset_index().rename(columns={'index':0,0:1})
		rnk_combined.columns = [0,1]
		rnk_combined = rnk_combined[rnk_combined[1]>0]
		print(rnk_combined)
		res = gp.prerank(rnk=rnk_combined,
						gene_sets='MSigDB_Hallmark_2020',
						processes=4,
						permutation_num=100,
						outdir=None,
						no_plot=True,seed=6).res2d.sort_values('fdr').iloc[0:5,3].reset_index()
		
		res = res[res['fdr']<0.01]
		
		if res.shape[0]>0:
			res['topic']=i
			df_grps = pd.concat([df_grps, res], axis=0, ignore_index=True)

	return df_grps

