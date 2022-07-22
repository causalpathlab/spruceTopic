import os
import pandas as pd
import math
import numpy as np
from torch import frac

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

def topic_top_genes(sp,top_n):

	top_genes = []
	top_genes = generate_gene_vals(sp.cell_topic.beta_mean,top_n,top_genes,'top_genes')

	return pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])

def topic_top_lr_genes(sp,top_n=5):
	top_genes = []
	top_genes = generate_gene_vals(sp.interaction_topic.beta_l,top_n,top_genes,'receptors')
	top_genes = generate_gene_vals(sp.interaction_topic.beta_r,top_n,top_genes,'ligands')

	return pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])
	

def topic_top_lr_pair_genes(sp,df_db,top_n=5):

	df_beta1 = sp.interaction_topic.beta_l
	df_beta2 = sp.interaction_topic.beta_r

	top_genes = []
	top_genes = generate_gene_vals(df_beta1,top_n,top_genes,'receptors')
	top_genes = generate_gene_vals(df_beta2,top_n,top_genes,'ligands')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])
	top_genes = df_top_genes['Gene'].unique()

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
	return df_lr_topic


def assign_gene_bias(args):
	args_home = os.environ['args_home']
	l_fname = args_home+args.input+args.raw_l_data_genes
	r_fname = args_home+args.input+args.raw_r_data_genes

	ligands = pd.read_pickle(l_fname)[0].values
	receptors = pd.read_pickle(r_fname)[0].values

	df_beta1_bias = pd.read_csv(args_home+args.output+args.interaction_model['out']+args.interaction_model['mfile']+'_ietm_beta1_bias.tsv.gz',sep='\t',compression='gzip')
	df_beta2_bias = pd.read_csv(args_home+args.output+args.interaction_model['out']+args.interaction_model['mfile']+'_ietm_beta2_bias.tsv.gz',sep='\t',compression='gzip')

	df_beta1_bias.columns = receptors
	df_beta2_bias.columns = ligands

	df_beta1_bias = df_beta1_bias.T.reset_index()
	df_beta1_bias.columns = ['gene','val']
	df_beta1_bias['group'] = 'mid'

	df_beta2_bias = df_beta2_bias.T.reset_index()
	df_beta2_bias.columns = ['gene','val']
	df_beta2_bias['group'] = 'mid'


	n = 25
	l25 = list(df_beta1_bias.sort_values('val').head(n)['gene'].values) 
	h25 = list(df_beta1_bias.sort_values('val',ascending=False).head(n)['gene'].values)
	df_beta1_bias.loc[df_beta1_bias['gene'].isin(l25),['group']] = 'l25'
	df_beta1_bias.loc[df_beta1_bias['gene'].isin(h25),['group']] = 'h25'
	l25 = list(df_beta2_bias.sort_values('val').head(n)['gene'].values) 
	h25 = list(df_beta2_bias.sort_values('val',ascending=False).head(n)['gene'].values)
	df_beta2_bias.loc[df_beta2_bias['gene'].isin(l25),['group']] = 'l25'
	df_beta2_bias.loc[df_beta2_bias['gene'].isin(h25),['group']] = 'h25'

	df_beta1_bias.to_csv(args_home+args.output+args.interaction_model['out']+args.interaction_model['mfile']+'_ietm_beta1_bias_v2.tsv.gz',sep='\t')
	df_beta2_bias.to_csv(args_home+args.output+args.interaction_model['out']+args.interaction_model['mfile']+'_ietm_beta2_bias_v2.tsv.gz',sep='\t')

def sample_cells_with_celltype(sp,cell_n=50):

	dfz = sp.cell_topic.h.copy()
	dfz.columns = [ x.replace('h','') for x in dfz.columns]

	dfz['label'] = [x.split('_')[len(x.split('_'))-1] for x in dfz['cell']]

	f='/home/BCCRC.CA/ssubedi/projects/data/GSE176078mix/GSE176078_metadata.csv.gz'
	dfl = pd.read_csv(f,compression='gzip')
	dfl = dfl.rename(columns={'Unnamed: 0':'cell'})

	dflabel = pd.DataFrame()
	dflabel['l1'] =  [x for x in dfz[dfz['label']=='GSE176078']['cell']]
	dflabel['l2'] =  [x.replace('_GSE176078','') for x in dfz[dfz['label']=='GSE176078']['cell']]
	dflabel = pd.merge(dflabel,dfl,right_on='cell',left_on='l2',how='left')
	label_index=8
	label = dfl.columns[label_index]
	dfz = pd.merge(dfz,dflabel[['l1',label]],right_on='l1',left_on='cell',how='left')
	dfz[label] = dfz[label].mask(dfz[label].isna(), dfz['label'])

	# df_h_sample = dfz.groupby(label).sample(frac=0.05, random_state=1)
	df_h_sample = dfz.groupby(label).sample(n=50, random_state=1)
	print(df_h_sample[label].value_counts())

	df_h_sample = df_h_sample.rename(columns={label:'Topic'})
	df_h_sample = df_h_sample.drop(columns=['label','l1'])
	return df_h_sample


def get_topics(spr,df_kmeans):
	df_h_state = spr.interaction_topic.neighbour_h.copy()
	df_h_state['interact_topic'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_h_state.iterrows()]
	dflatent = pd.merge(df_h_state[['cell','interact_topic']],df_kmeans,how='left',on='cell')
	return dflatent

def add_interaction_topics(spr,df_kmeans):
	df_h_state = spr.interaction_topic.neighbour_h.copy()
	df_h_state['it_topic'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_h_state.iterrows()]

	df_kmeans = pd.merge(df_kmeans,df_h_state[['cell','it_topic']],how='left',on='cell')
	
	return df_kmeans

def topics_summary(spr,df_kmeans):

	df_h_celltype = spr.cell_topic.h.copy()
	df_h_celltype['cell_topic'] = df_h_celltype.iloc[:,1:].idxmax(axis=1)

	df_h_state = spr.interaction_topic.neighbour_h.copy()
	df_h_state['interact_topic'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_h_state.iterrows()]

	dflatent = pd.merge(df_h_state[['cell','interact_topic']],df_h_celltype[['cell','cell_topic']],how='left',on='cell')

	dflatent = pd.merge(dflatent,df_kmeans[['cell','cluster']],how='left',on='cell')

	dfsummary = dflatent.groupby(['cluster','cell_topic','interact_topic']).count()
	dfsummary = dfsummary.reset_index()
	
	return dfsummary


def get_cell_neighbours_states(df_nbr,df_its,cells):

	cc_idxs = df_nbr[df_nbr['cell'].isin(cells)].index

	res = []
	for idx in cc_idxs:
		cell = df_nbr.iloc[idx,0]
		nbr_idxs = df_nbr.iloc[idx,1:].values
		nbr_states = df_its.iloc[idx,1:].values
		nbr_cells = df_its.iloc[nbr_idxs]['cell'].values
		res.append([cell,nbr_cells,nbr_states])

	return res

def plt_cn(spr,df,statelist):

	import matplotlib.pylab as plt
	import colorcet as cc
	import seaborn as sns
	plt.rcParams['figure.figsize'] = [16,12]
	plt.rcParams['figure.autolayout'] = True

	fig, ax = plt.subplots(7,2)
 
	for idx,state in enumerate(statelist):
		df_sel = df[df['interact_topic'] == state]
		dfl,dfr = get_cell_neighbours_states_lr(spr,df_sel,state)
		sns.heatmap(dfl,annot=False,ax=ax[idx,0],cmap='viridis',yticklabels=['c','n'])
		ax[idx,0].set_ylabel(state)
		# ax[0].set_title('Ligands_interaction_topic_'+str(state))
		sns.heatmap(dfr,annot=False,ax=ax[idx,1],cmap='viridis',yticklabels=['c','n'])
		# ax[1].set_title('Receptors_interaction_topic_'+str(state))
	plt.savefig(spr.interaction_topic.model_id+'_lr_interaction_topic_hmap.pdf')
	plt.close()

def get_cell_neighbours_states_lr(spr,df,state):
	df_cancer_l = spr.data.raw_l_data[spr.data.raw_l_data['index'].isin(df['Cancer'].values)]
	df_cancer_r = spr.data.raw_r_data[spr.data.raw_r_data['index'].isin(df['Cancer'].values)]

	df_nbr_l = spr.data.raw_l_data[spr.data.raw_l_data['index'].isin(df['nbr'].values)]
	df_nbr_r = spr.data.raw_r_data[spr.data.raw_r_data['index'].isin(df['nbr'].values)]

	# df_cancer_l.iloc[:,1:] = df_cancer_l.iloc[:,1:].div(df_cancer_l.iloc[:,1:].sum(axis=1), axis=0)
	# df_cancer_r.iloc[:,1:] = df_cancer_r.iloc[:,1:].div(df_cancer_r.iloc[:,1:].sum(axis=1), axis=0)

	# df_nbr_l.iloc[:,1:] = df_nbr_l.iloc[:,1:].div(df_nbr_l.iloc[:,1:].sum(axis=1), axis=0)
	# df_nbr_r.iloc[:,1:] = df_nbr_r.iloc[:,1:].div(df_nbr_r.iloc[:,1:].sum(axis=1), axis=0)

	# df_cancer_l = df_cancer_l.sample(n=df_nbr_l.shape[0])
	# df_cancer_r = df_cancer_r.sample(n=df_nbr_r.shape[0])


	top_r = list(spr.interaction_topic.beta_r.iloc[state,:].sort_values(ascending=False).head(25).index)
	top_l = list(spr.interaction_topic.beta_l.iloc[state,:].sort_values(ascending=False).head(25).index)
	
	df_cancer_l = df_cancer_l.loc[:,top_l]
	df_cancer_r = df_cancer_r.loc[:,top_r]
	df_nbr_l = df_nbr_l.loc[:,top_l]
	df_nbr_r = df_nbr_r.loc[:,top_r]

	df_cancer_l = df_cancer_l.mean()
	df_cancer_r = df_cancer_r.mean()
	df_nbr_l = df_nbr_l.mean()
	df_nbr_r = df_nbr_r.mean()

	dfl = pd.DataFrame([df_cancer_l,df_nbr_l],index=['cancer','nbr'])
	dfr = pd.DataFrame([df_cancer_r,df_nbr_r],index=['cancer','nbr'])

	dfl += 1e-5
	dfr += 1e-5
	dfl = dfl.applymap(math.sqrt)
	dfr = dfr.applymap(math.sqrt)
	
	return dfl,dfr

def interaction_summary(sp,topics_prob):
	summary = []
	for idx in range(len(topics_prob)):
		ci = sp.data.raw_l_data['index'][idx]
		for nbr_idx in range(len(topics_prob[idx])):
			cj = sp.data.raw_l_data['index'][sp.data.neighbour.iloc[idx,nbr_idx]]
			for topic_idx,interaction_p in enumerate(topics_prob[idx][nbr_idx]):
				summary.append([ci,cj,topic_idx,interaction_p])
		if idx % 10000 == 0:
			print(idx)
	df = pd.DataFrame(summary,columns=['cell_i','cell_j','interaction_topic','interaction_probability'])
	df.to_csv(sp.interaction_topic.model_id+'_interaction_summary.tsv.gz',sep='\t',index=False,compression='gzip')
