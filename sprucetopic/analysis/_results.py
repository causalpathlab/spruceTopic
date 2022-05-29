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

def topic_top_genes(args,top_n=5):

	args_home = os.environ['args_home']

	df_genes=pd.read_pickle(args_home+args.input+args.raw_data_genes)

	df_beta = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_beta.tsv.gz',sep='\t',compression='gzip')

	df_beta.columns = df_genes[0].values

	top_genes = []
	top_genes = generate_gene_vals(df_beta,top_n,top_genes,'top_genes')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])

	df_top_genes.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_top_'+str(top_n)+'_genes_topic.tsv.gz',sep='\t',index=False)

def topic_top_lr_genes(args,top_n=5):

	args_home = os.environ['args_home']
	l_fname = args_home+args.input+args.raw_l_data_genes
	r_fname = args_home+args.input+args.raw_r_data_genes

	ligands = pd.read_pickle(l_fname)[0].values
	receptors = pd.read_pickle(r_fname)[0].values

	df_beta1 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_beta1.tsv.gz',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_beta2.tsv.gz',sep='\t',compression='gzip')

	df_beta1.columns = receptors
	df_beta2.columns = ligands

	top_genes = []
	top_genes = generate_gene_vals(df_beta1,top_n,top_genes,'receptors')
	top_genes = generate_gene_vals(df_beta2,top_n,top_genes,'ligands')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])
	df_top_genes.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_top_'+str(top_n)+'_genes_topic.tsv.gz',sep='\t',index=False)

def topic_top_lr_pair_genes(args,top_n=5):

	args_home = os.environ['args_home']
	l_fname = args_home+args.input+args.raw_l_data_genes
	r_fname = args_home+args.input+args.raw_r_data_genes

	ligands = pd.read_pickle(l_fname)[0].values
	receptors = pd.read_pickle(r_fname)[0].values

	df_beta1 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_beta1.tsv.gz',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_beta2.tsv.gz',sep='\t',compression='gzip')

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
	df_lr_topic.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'top_'+str(top_n)+'_lrpair_topic.tsv.gz',sep='\t')

def assign_gene_bias(args):
	args_home = os.environ['args_home']
	l_fname = args_home+args.input+args.raw_l_data_genes
	r_fname = args_home+args.input+args.raw_r_data_genes

	ligands = pd.read_pickle(l_fname)[0].values
	receptors = pd.read_pickle(r_fname)[0].values

	df_beta1_bias = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_beta1_bias.tsv.gz',sep='\t',compression='gzip')
	df_beta2_bias = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_beta2_bias.tsv.gz',sep='\t',compression='gzip')

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



	df_beta1_bias.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_beta1_bias_v2.tsv.gz',sep='\t')
	df_beta2_bias.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_beta2_bias_v2.tsv.gz',sep='\t')

	# pd.DataFrame(pd.cut(df_beta1_bias[0],bins=10).reset_index().groupby(0)['index'].apply(','.join)).to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta1_bias_receptors_table.tsv.gz',sep='\t')
	# pd.DataFrame(pd.cut(df_beta2_bias[0],bins=10).reset_index().groupby(0)['index'].apply(','.join)).to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta2_bias_ligands_table.tsv.gz',sep='\t')

def sample_cells_with_latent(args,cell_n=50):

	args_home = os.environ['args_home']
	dfh = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_h.tsv.gz',sep='\t',compression='gzip')
	dfh.columns = [ x.replace('h','') for x in dfh.columns]
	# dfh['Topic'] = dfh.iloc[:,1:].idxmax(axis=1)

	# df_h_sample = dfh.groupby('Topic').sample(frac=0.1, random_state=1)
	# df_h_sample.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_h_topic_sample.tsv',sep='\t',index=False)

	dfz=dfh
	dfz['label'] = [x.split('_')[len(x.split('_'))-1] for x in dfz['cell']]

	f='/home/sishirsubedi/projects/spruce_topic/input/GSEmix/GSE176078_metadata.csv.gz'
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
	df_h_sample.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_h_topic_sample.tsv',sep='\t',index=False)

def topics_summary(args):
	
	args_home = os.environ['args_home']
	model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']

	df_h_celltype = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_h.tsv.gz',sep='\t',compression='gzip')

	df_h_state = pd.read_csv(model_file+'_ietm_interaction_states.tsv.gz',sep='\t',compression='gzip')

	df_h_state['state'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_h_state.iterrows()]

	df_h_celltype['topic'] = df_h_celltype.iloc[:,1:].idxmax(axis=1)

	dflatent = pd.merge(df_h_state[['cell','state']],df_h_celltype[['cell','topic']],how='left',on='cell')

	dflatent['label'] = [x.split('_')[len(x.split('_'))-1] for x in dflatent['cell']]

	f=args_home+args.input+'GSE176078_metadata.csv.gz'
	dfm = pd.read_csv(f,compression='gzip')
	dfm = dfm.rename(columns={'Unnamed: 0':'gse_cell'})

	dflabel = pd.DataFrame()
	dflabel['cell'] =  [x for x in dflatent[dflatent['label']=='GSE176078']['cell']]
	dflabel['gse_cell'] =  [x.replace('_GSE176078','') for x in dflatent[dflatent['label']=='GSE176078']['cell']]
	dflabel = pd.merge(dflabel,dfm,on='gse_cell',how='left')

	dflatent = pd.merge(dflatent,dflabel,on='cell',how='left')
	dflatent = dflatent[dflatent['label']=='GSE176078']

	dflabel.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_cell_label_topic_state_meta.tsv.gz',sep='\t',index=False,compression='gzip')

def topics_summary(args):
	
	args_home = os.environ['args_home']
	model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']

	df_h_celltype = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_netm_h.tsv.gz',sep='\t',compression='gzip')

	df_h_state = pd.read_csv(model_file+'_ietm_interaction_states.tsv.gz',sep='\t',compression='gzip')

	df_h_state['state'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_h_state.iterrows()]

	df_h_celltype['topic'] = df_h_celltype.iloc[:,1:].idxmax(axis=1)

	dflatent = pd.merge(df_h_state[['cell','state']],df_h_celltype[['cell','topic']],how='left',on='cell')

	dflatent['label'] = [x.split('_')[len(x.split('_'))-1] for x in dflatent['cell']]

	# f=args_home+args.input+'GSE176078_metadata.csv.gz'
	f='/home/sishirsubedi/projects/spruce_topic/input/GSEmix/GSE176078_metadata.csv.gz'
	dfl = pd.read_csv(f,compression='gzip')
	dfl = dfl.rename(columns={'Unnamed: 0':'cell'})

	dflabel = pd.DataFrame()
	dflabel['l1'] =  [x for x in dflatent[dflatent['label']=='GSE176078']['cell']]
	dflabel['l2'] =  [x.replace('_GSE176078','') for x in dflatent[dflatent['label']=='GSE176078']['cell']]
	dflabel = pd.merge(dflabel,dfl,right_on='cell',left_on='l2',how='left')

	label_index = 8
	label = dfl.columns[label_index]

	dflatent = pd.merge(dflatent,dflabel[['l1',label]],right_on='l1',left_on='cell',how='left')
	dflatent[label] = dflatent[label].mask(dflatent[label].isna(), dflatent['label'])

	dfsummary = dflatent.groupby([label,'topic','state']).count()
	dfsummary = dfsummary.reset_index()
	dfsummary = dfsummary.iloc[:,:-2]

	dfsummary = dfsummary.rename(columns={label:'label'})

	dfsummary.to_csv(model_file+'model_summary.csv',index=False)

def cell_interaction_network(args):
	
	import igraph
	import seaborn as sns

	args_home = os.environ['args_home']

	model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']

	df = pd.read_csv(model_file+'_ietm_interaction_states.tsv.gz',sep='\t',compression='gzip')
	df['state'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df.iterrows()]
	df['label'] = [x.split('_')[len(x.split('_'))-1] for x in df['cell']]

	# meta_path = args_home+ args.input+args.metadata
	f='/home/sishirsubedi/projects/spruce_topic/input/GSEmix/GSE176078_metadata.csv.gz'
	df_meta = pd.read_csv(f)
	df_meta = df_meta.rename(columns={'Unnamed: 0':'cell'})

	dflabel = pd.DataFrame()
	dflabel['l1'] =  [x for x in df[df['label']=='GSE176078']['cell']]
	dflabel['l2'] =  [x.replace('_GSE176078','') for x in df[df['label']=='GSE176078']['cell']]

	print(dflabel.columns)
	print(df_meta.columns)
	dflabel = pd.merge(dflabel,df_meta,right_on='cell',left_on='l2',how='left')
	label = df_meta.columns[8]

	df = pd.merge(df,dflabel[['l1',label]],right_on='l1',left_on='cell',how='left')
	df[label] = df[label].mask(df[label].isna(), df['label'])

	df = df.dropna()
	dfs = df.groupby('label', group_keys=False).apply(pd.DataFrame.sample,frac=0.1)
	# dfs = df.groupby('label', group_keys=False).apply(pd.DataFrame.sample,100)

	dfs = dfs.drop(columns=['l1','label'])
	dfs = dfs.rename(columns={label:'label'})

	cells = [x for x in list(dfs.label.values)]

	states = ['state'+str(int(x)) for x in dfs['state'].unique()]
	print(states)

	## construct graph
	# if len(states)<=10:
	# 	colors = sns.color_palette("Paired")
	# elif len(states)>10:
	colors = sns.color_palette('Spectral',int(dfs['state'].unique().max())+1)

	g = igraph.Graph()
	g.add_vertices(cells+states)

	gedges = []
	for i,row in dfs.iterrows():
		cell = row.label
		s = int(row['state'])
		interaction_state = 'state'+str(s)
		epair = cell + '_' + interaction_state

		if epair not in gedges:
			g.add_edge(cell,interaction_state,weight=0.1,color=colors[int(s)])
			gedges.append(epair)
		else:
			eid = g.get_eid(cell,interaction_state)
			g.es[eid]['weight'] += 0.1


	to_delete_eids = [e.index for e in g.es if e['weight']<5]
	g.delete_edges(to_delete_eids)

	to_delete_vids = [v.index for v in g.vs if v.degree()<1]
	g.delete_vertices(to_delete_vids)

	
	for v in g.vs:
		if "state" not in v['name']: 
			v['attr'] = 'cell'
			v['size'] = 50
			v['color'] = 'deeppink2'
			v['shape'] = 'circle'
		else: 
			v['attr'] = 'state'
			v['size'] = 20
			v['color'] = colors[int(str(v['name']).replace('state',''))]
			v['shape'] = 'square'
		

	visual_style={}
	visual_style["vertex_label"] = [ x.replace('state','s') for x in g.vs["name"]]
	visual_style["vertex_size"] = g.vs['size']
	visual_style["vertex_color"] = g.vs['color']
	visual_style["vertex_label_size"] = 12
	visual_style["edge_color"] = g.es['color']
	visual_style["edge_width"] = [x/7 for x in g.es['weight']]
	visual_style["layout"] = g.layout_kamada_kawai()
	igraph.plot(g, target=model_file+"_network_10.pdf",**visual_style,margin=40)

def cancer_centric_view_lr(args):

	args_home = os.environ['args_home']
	l_fname = args_home+args.input+args.raw_l_data_genes
	r_fname = args_home+args.input+args.raw_r_data_genes

	ligands = pd.read_pickle(l_fname)[0].values
	receptors = pd.read_pickle(r_fname)[0].values

	df_beta1 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_beta1.tsv.gz',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_beta2.tsv.gz',sep='\t',compression='gzip')

	df_beta1.columns = receptors
	df_beta2.columns = ligands

	top_r_genes = []
	top_r_genes = generate_gene_vals(df_beta1.iloc[[5],:],10,top_r_genes,'receptors')
	top_r_genes = pd.Series([x[3] for x in top_r_genes]).unique()

	top_l_genes = []
	top_l_genes = generate_gene_vals(df_beta2.iloc[[5],:],10,top_l_genes,'ligands')
	top_l_genes = pd.Series([x[3] for x in top_l_genes]).unique()


	dfm = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_cell_label_topic_state_meta.tsv.gz',sep='\t',compression='gzip')
	# cancer_cells = dfm[dfm.celltype_major=='Cancer Epithelial']['cell'].values
	cancer_cells = dfm[dfm.celltype_major.isin(['Cancer Epithelial','T-cells','B-cells'])]['cell'].values

	df_h_state = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'_ietm_interaction_states.tsv.gz',sep='\t',compression='gzip')
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

	model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']
	df_l.to_csv(model_file+'_s5_cv_ligands.tsv.gz',index=False,sep='\t',compression='gzip')
	df_r.to_csv(model_file+'_s5_cv_receptors.tsv.gz',index=False,sep='\t',compression='gzip')
