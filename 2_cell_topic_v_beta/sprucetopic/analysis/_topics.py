import os
import pandas as pd
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

	sp.interaction_topic.beta1.columns = sp.data.raw_r_data_genes
	sp.interaction_topic.beta2.columns = sp.data.raw_l_data_genes

	top_genes = []
	top_genes = generate_gene_vals(sp.interaction_topic.beta1,top_n,top_genes,'receptors')
	top_genes = generate_gene_vals(sp.interaction_topic.beta2,top_n,top_genes,'ligands')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])
	df_top_genes.to_csv(sp.interaction_topic.model_id+'_ietm_top_'+str(top_n)+'_genes_topic.tsv.gz',sep='\t',index=False)

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
	df_h_state = spr.interaction_topic.neighbour_h
	df_h_state['state'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_h_state.iterrows()]
	dflatent = pd.merge(df_h_state[['cell','state']],df_kmeans,how='left',on='cell')
	return dflatent


def topics_summary(spr,df_kmeans):

	df_h_celltype = spr.cell_topic.h
	df_h_celltype['topic'] = df_h_celltype.iloc[:,1:].idxmax(axis=1)

	df_h_state = spr.interaction_topic.neighbour_h
	df_h_state['state'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_h_state.iterrows()]

	dflatent = pd.merge(df_h_state[['cell','state']],df_h_celltype[['cell','topic']],how='left',on='cell')

	dflatent = pd.merge(dflatent,df_kmeans[['cell','cluster']],how='left',on='cell')

	dfsummary = dflatent.groupby(['cluster','topic','state']).count()
	dfsummary = dfsummary.reset_index()
	
	return dfsummary


def cancer_centric_view_lr(args):

	args_home = os.environ['args_home']
	l_fname = args_home+args.input+args.raw_l_data_genes
	r_fname = args_home+args.input+args.raw_r_data_genes

	ligands = pd.read_pickle(l_fname)[0].values
	receptors = pd.read_pickle(r_fname)[0].values

	df_beta1 = pd.read_csv(args_home+args.output+args.interaction_model['out']+args.interaction_model['mfile']+'_ietm_beta1.tsv.gz',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.interaction_model['out']+args.interaction_model['mfile']+'_ietm_beta2.tsv.gz',sep='\t',compression='gzip')

	df_beta1.columns = receptors
	df_beta2.columns = ligands

	top_r_genes = []
	top_r_genes = generate_gene_vals(df_beta1.iloc[[5],:],10,top_r_genes,'receptors')
	top_r_genes = pd.Series([x[3] for x in top_r_genes]).unique()

	top_l_genes = []
	top_l_genes = generate_gene_vals(df_beta2.iloc[[5],:],10,top_l_genes,'ligands')
	top_l_genes = pd.Series([x[3] for x in top_l_genes]).unique()


	dfm = pd.read_csv(args_home+args.output+args.interaction_model['out']+args.interaction_model['mfile']+'_ietm_cell_label_topic_state_meta.tsv.gz',sep='\t',compression='gzip')
	# cancer_cells = dfm[dfm.celltype_major=='Cancer Epithelial']['cell'].values
	cancer_cells = dfm[dfm.celltype_major.isin(['Cancer Epithelial','T-cells','B-cells'])]['cell'].values

	df_h_state = pd.read_csv(args_home+args.output+args.interaction_model['out']+args.interaction_model['mfile']+'_ietm_interaction_states.tsv.gz',sep='\t',compression='gzip')
	df_h_state['state'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_h_state.iterrows()]
	df_h_state = df_h_state[['cell','state']]
	# df_h_state = df_h_state[df_h_state['cell'].isin(cancer_cells)]

	l_fname = args_home+args.input+args.raw_l_data
	r_fname = args_home+args.input+args.raw_r_data

	df_l = pd.read_pickle(l_fname)
	df_r = pd.read_pickle(r_fname)

	df_l = df_l[df_l['index'].isin(cancer_cells)]
	df_r = df_r[df_r['index'].isin(cancer_cells)]

	df_l = pd.merge(df_l,df_h_state,right_on='cell',left_on='index',how='left')
	df_r = pd.merge(df_r,df_h_state,right_on='cell',left_on='index',how='left')


	df_l = df_l[['cell','state']+list(top_l_genes)]
	df_r = df_r[['cell','state']+list(top_r_genes)]

	model_file = args_home+args.output+args.interaction_model['out']+args.interaction_model['mfile']
	df_l.to_csv(model_file+'_s5_cv_ligands.tsv.gz',index=False,sep='\t',compression='gzip')
	df_r.to_csv(model_file+'_s5_cv_receptors.tsv.gz',index=False,sep='\t',compression='gzip')

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
