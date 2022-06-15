from audioop import add
import os
from turtle import title
import pandas as pd
import numpy as np
import igraph
import seaborn as sns


def cell_interaction_network(args):
	

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

def lr_network(sp):
	for i in range(25):	
		df = pd.read_csv(sp.model_id+'topic_'+str(i)+'_ietm_alpha.tsv.gz',sep='\t',compression='gzip')
		
		# df.index = sp.data.raw_l_data_genes
		# df.columns = sp.data.raw_r_data_genes

		df.columns = sp.data.raw_l_data_genes
		df.index = sp.data.raw_r_data_genes

		th = 0.99
		g = igraph.Graph()
		# g.add_vertices(sp.data.raw_l_data_genes + sp.data.raw_r_data_genes)

		added = []
		for ligand in df.index:
			for receptor in df.columns:
				w = df.loc[ligand,receptor]
				if w >= th:
					if ligand not in added:
						g.add_vertex(ligand)
						added.append(ligand)
					if receptor not in added:
						g.add_vertex(receptor)
						added.append(receptor)
					g.add_edge(ligand,receptor,weight=w)
		
		# to_delete_vids = [v.index for v in g.vs if v.degree()<1]
		# g.delete_vertices(to_delete_vids)

		for v in g.vs:
			if v['name'] in list(df.columns): 
				v['attr'] = 'receptor'
				v['size'] = 50
				v['color'] = 'deeppink2'
				v['shape'] = 'circle'
			else: 
				v['attr'] = 'ligand'
				v['size'] = 20
				v['color'] = 'skyblue'
				v['shape'] = 'square'
			

		visual_style={}
		visual_style["vertex_label"] =  g.vs['name']
		visual_style["vertex_size"] = g.vs['size']
		visual_style["vertex_color"] = g.vs['color']
		visual_style["vertex_label_size"] = 12
		visual_style["layout"] = g.layout_kamada_kawai()
		igraph.plot(g, target=sp.model_id+'_lr_network_topic'+str(i)+'.pdf',**visual_style,margin=40,title='topic_'+str(i))
