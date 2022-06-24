from audioop import add
import os
from turtle import title
import pandas as pd
import numpy as np
import igraph
import seaborn as sns


def cell_interaction_network(spr,df):

	dfs = df.groupby('cluster', group_keys=False).apply(pd.DataFrame.sample,frac=0.1)
	dfs['cluster'] = [str(x)+'-'+y for x,y in zip(dfs['cluster'],dfs['cluster_celltype'])]


	cells = [x for x in list(dfs['cluster'].values)]

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
		cell = row.cluster
		s = int(row['state'])
		interaction_state = 'state'+str(s)
		epair = cell + '_' + interaction_state

		if epair not in gedges:
			g.add_edge(cell,interaction_state,weight=0.1,color=colors[int(s)])
			gedges.append(epair)
		else:
			eid = g.get_eid(cell,interaction_state)
			g.es[eid]['weight'] += 0.1


	to_delete_eids = [e.index for e in g.es if e['weight']<10]
	g.delete_edges(to_delete_eids)

	to_delete_vids = [v.index for v in g.vs if v.degree()<1]
	g.delete_vertices(to_delete_vids)

	
	for v in g.vs:
		if "state" not in v['name']: 
			v['attr'] = 'cell'
			v['size'] = 10
			v['color'] = 'aqua'
			v['shape'] = 'circle'
		else: 
			v['attr'] = 'state'
			v['size'] = 20
			v['color'] = colors[int(str(v['name']).replace('state',''))]
			v['shape'] = 'square'
		

	visual_style={}
	visual_style["vertex_label"] = [ x.replace('state','s') for x in g.vs["name"]]
	# visual_style["vertex_size"] = g.vs['size']
	# visual_style["vertex_color"] = g.vs['color']
	visual_style["vertex_label_size"] = 10
	visual_style["edge_color"] = g.es['color']
	visual_style["edge_width"] = [x/7 for x in g.es['weight']]
	# visual_style["layout"] = g.layout_circle()
	visual_style["layout"] = g.layout_reingold_tilford_circular()
	igraph.plot(g, target=spr.interaction_topic.model_id+"_it_model_network_plot.png",**visual_style,margin=40)

def interaction_statewise_lr_network(spr,states,df_db=None):

	import matplotlib.pylab as plt
	import colorcet as cc
	import seaborn as sns
	plt.rcParams['figure.figsize'] = [20,10]
	plt.rcParams['figure.autolayout'] = True

	fig, ax = plt.subplots(2,4) 
	ax = ax.ravel()

	ligands = list(spr.interaction_topic.beta_r.columns) 
	receptors = list(spr.interaction_topic.beta_l.columns) 

	for i,state in enumerate(states): 	
		l_vals =  spr.interaction_topic.beta_r.iloc[state,:].values
		r_vals =  spr.interaction_topic.beta_l.iloc[state,:].values

		l = l_vals.reshape(l_vals.shape[0],1)
		r = r_vals.reshape(r_vals.shape[0],1)
		r = r.T
		lr = l[:,None]+r
		lr = lr.flatten()
		lr = lr[lr.argsort()]
		cuttoff = lr[-25:][0]
		print('cuttoff is...'+str(cuttoff))


		# lr_topic = []
		# lr_pair = []
		# for lr_p in df_db['lr_pair']:
		# 	l,r = lr_p.split('_')[0],lr_p.split('_')[1]
		# 	if l in df_beta2.columns and r in df_beta1.columns and \
		# 		l in top_genes and r in top_genes:
		# 		lr_topic.append((df_beta2[l]+df_beta1[r])/2)
		# 		lr_pair.append(lr_p)
		# df_lr_topic = pd.DataFrame(lr_topic)
		# df_lr_topic.index=lr_pair
		
		g = igraph.Graph()

		weights = []
		for li,ligand in enumerate(ligands):
			for ri,receptor in enumerate(receptors):
				w = l_vals[li]+r_vals[ri]
				if w>cuttoff:
					g.add_vertex(ligand)
					g.add_vertex(receptor)
					g.add_edge(ligand,receptor,weight= w)
		

		# to_delete_eids = [e.index for e in g.es if e['weight']<10]
		# g.delete_edges(to_delete_eids)

		to_delete_vids = [v.index for v in g.vs if v.degree()<1]
		g.delete_vertices(to_delete_vids)

		for v in g.vs:
			if v['name'] in receptors: 
				v['attr'] = 'receptor'
				v['size'] = 15
				v['color'] = 'r'
				v['shape'] = 'diamond'
			else: 
				v['attr'] = 'ligand'
				v['size'] = 10
				v['color'] = 'b'
				v['shape'] = 'square'
			

		visual_style={}
		visual_style["vertex_label"] =  g.vs['name']
		visual_style["vertex_size"] = g.vs['size']
		visual_style["vertex_color"] = g.vs['color']
		visual_style["vertex_label_size"] = 14
		# visual_style["edge_width"] = [x/5 for x in g.es['weight']]
		visual_style["layout"] = g.layout_circle()

		ax[i].set_axis_off()
		ax[i].set_title('interaction topic_'+str(state))
		igraph.plot(g,target=ax[i],**visual_style,margin=20,title='interaction topic_'+str(state))
		# igraph.plot(g, target=spr.interaction_topic.model_id+'_lr_network_it_state_'+str(state)+'.png',**visual_style,margin=40,)
	plt.savefig(spr.interaction_topic.model_id+'_lr_network_it_states.png');plt.close()
