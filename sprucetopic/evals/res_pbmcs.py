import os
import pandas as pd
import numpy as np

cell_type={
		0:'CD4T-Naive',
		1:'CD4T-Memory',
		4:'CD8T',
		6:'NK',
		3:'B',
		2:'CD14',
		5:'FCGR3A',
		7:'Dendritic',
		8:'Megakaryocytes'
	}

def pbmc_sample_cells_with_latent(args,cell_n=50):

	args_home = os.environ['args_home']

	df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_hh_data.tsv',sep='\t',compression='gzip')

	meta_path = args_home+args.scanpy_output+'pbmc_cluster.pkl'
	df_meta = pd.read_pickle(meta_path)
	dfjoin = pd.merge(df_h,df_meta,right_on='index',left_on='cell',how='left')
	dfjoin = dfjoin[dfjoin['leiden'].notna()]


	# df_h_sample = dfjoin.groupby('leiden').sample(cell_n, random_state=1)
	df_h_sample = dfjoin.drop(columns='index').rename(columns={'leiden':'cluster'})
	df_h_sample['cluster'] = [cell_type[int(x)] for x in df_h_sample['cluster']]
	
	df_h_sample.columns = [ x.replace('hh','k') for x in df_h_sample.columns]
	df_h_sample.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'hh_cell_topic_sample.tsv',sep='\t',index=False)



def topic_celltype_marker_genes(args):

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

	marker_genes = ['IL7R', 'S100A4','CCR7','CD79A', 'MS4A1', 'CD8A', 'CD8B', 
				   'LYZ', 'CD14',
				'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
				'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
	
	# In [43]: [x for x in ligand_receptors if x in marker_genes]
	# Out[43]: ['IL7R', 'FCER1A', 'CD79A', 'KLRB1', 'LYZ', 'PPBP', 'CD14', 'S100A8']
	
	df_beta = pd.concat([df_beta1,df_beta2],axis=1)

	marker_genes = [x for x in marker_genes if x in df_beta.columns]

	df_beta = df_beta[marker_genes]

	df_beta.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'marker_genes_topic.tsv',sep='\t',index=False)

def topic_celltype_marker_genes_lrnet(args):

	args_home = os.environ['args_home']

	df_l_genes=pd.read_pickle(args_home+args.input+args.raw_l_data_genes)
	df_r_genes=pd.read_pickle(args_home+args.input+args.raw_r_data_genes)


	df_beta1 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta1_data.tsv',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta2_data.tsv',sep='\t',compression='gzip')

	df_beta1.columns = df_r_genes[0].values
	df_beta2.columns = df_l_genes[0].values

	marker_genes = ['IL7R', 'S100A4','CCR7','CD79A', 'MS4A1', 'CD8A', 'CD8B', 
				   'LYZ', 'CD14',
				'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
				'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
	
	# In [43]: [x for x in ligand_receptors if x in marker_genes]
	# Out[43]: ['IL7R', 'FCER1A', 'CD79A', 'KLRB1', 'LYZ', 'PPBP', 'CD14', 'S100A8']
	
	df_beta = pd.concat([df_beta1,df_beta2],axis=1)

	marker_genes = [x for x in marker_genes if x in df_beta.columns]

	df_beta = df_beta[marker_genes]

	df_beta.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'marker_genes_topic.tsv',sep='\t',index=False)

def combine_topics_pbmc(args):
	
	args_home = os.environ['args_home']
	model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']

	df_h_celltype = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_hh_data.tsv',sep='\t',compression='gzip')

	df_h_state = pd.read_csv(model_file+'lrnet_interaction_states.tsv',sep='\t',compression='gzip')

	df_h_state['state'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_h_state.iterrows()]

	df_h_celltype['topic'] = df_h_celltype.iloc[:,1:].idxmax(axis=1)

	dflatent = pd.merge(df_h_state[['cell','state']],df_h_celltype[['cell','topic']],how='left',on='cell')

	meta_path = args_home+args.scanpy_output+'pbmc_cluster.pkl'
	df_meta = pd.read_pickle(meta_path)

	dflatent['label'] = pd.merge(dflatent,df_meta,right_on='index',left_on='cell',how='left')['leiden'].values

	dfsummary = dflatent.groupby(['label','topic','state']).count()
	dfsummary = dfsummary.reset_index()
	dfsummary['label'] = [cell_type[int(x)] for x in dfsummary['label']]
	dfsummary.to_csv(model_file+'model_summary.csv',index=False)

def cell_interaction_network(args):
	
	import igraph 
	import seaborn as sns

	args_home = os.environ['args_home']

	model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']

	df = pd.read_csv(model_file+'lrnet_interaction_states.tsv',sep='\t',compression='gzip')

	# df = df.sample(frac=1).reset_index(drop=True)

	meta_path = args_home+args.scanpy_output+'pbmc_cluster.pkl'
	df_meta = pd.read_pickle(meta_path)
	df['label'] = pd.merge(df,df_meta,right_on='index',left_on='cell',how='left')['leiden'].values
	df = df.dropna()
	df['label'] = [cell_type[int(x)] for x in df['label']]

	df['state'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df.iterrows()]

	dfs = df.groupby('label', group_keys=False).apply(pd.DataFrame.sample, frac=1)
	cells = [x for x in list(dfs.label.values)]

	states = ['t'+str(x) for x in dfs['state'].unique()]
	print(states)
	## construct graph
	if len(states)<=10:
		colors = sns.color_palette("Paired")
	elif len(states)>10:
		colors = sns.color_palette('Spectral',dfs['state'].unique().max()+1)
	print(colors)

	g = igraph.Graph()
	g.add_vertices(cells+states)

	gedges = []
	for i,row in dfs.iterrows():
		cell = row.label
		s = row['state']
		interaction_state = 't'+str(s)
		epair = cell + '_' + interaction_state

		if epair not in gedges:
			g.add_edge(cell,interaction_state,weight=0.1,color=colors[s])
			gedges.append(epair)
		else:
			eid = g.get_eid(cell,interaction_state)
			g.es[eid]['weight'] += 0.1


	to_delete_eids = [e.index for e in g.es if e['weight']<1 ]
	g.delete_edges(to_delete_eids)

	to_delete_vids = [v.index for v in g.vs if v.degree()<1]
	g.delete_vertices(to_delete_vids)

	
	for v in g.vs:
		if v['name'] in cell_type.values(): 
			v['attr'] = 'cell'
			v['size'] = 50
			v['color'] = 'deeppink2'
			v['shape'] = 'circle'

		else: 
			v['attr'] = 'state'
			v['size'] = 20
			v['color'] = colors[int(v['name'][1])]
			v['shape'] = 'square'
		

	visual_style={}
	visual_style["vertex_label"] = g.vs["name"]
	visual_style["vertex_size"] = g.vs['size']
	visual_style["vertex_color"] = g.vs['color']
	visual_style["vertex_label_size"] = 12
	visual_style["edge_color"] = g.es['color']
	visual_style["edge_width"] = [x/7 for x in g.es['weight']]
	visual_style["layout"] = g.layout_kamada_kawai()
	igraph.plot(g, target=model_file+"network.pdf",**visual_style,margin=40)




