import os
import pandas as pd
import numpy as np
from premodel import prep_tcells as prep

def topic_top_genes(args,top_n=5):

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

	top_genes = []
	top_genes = prep.generate_gene_vals(df_beta1,top_n,top_genes,'T cell')
	top_genes = prep.generate_gene_vals(df_beta2,top_n,top_genes,'non T cell')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])

	db = args_home+args.database+args.tcell_signature_genes

	df_cd4 = pd.read_excel(db,sheet_name='CD4',header=1,usecols='A,B,C' )
	df_cd4 = df_cd4.sort_values('comb.ES',ascending=False).drop_duplicates('geneSymbol')

	df_cd8 = pd.read_excel(db,sheet_name='CD8',header=1,usecols='A,B,C' )
	df_cd8 = df_cd8.sort_values('comb.ES',ascending=False).drop_duplicates('geneSymbol')

	df_cd48_top_genes = pd.DataFrame(df_top_genes['Gene'].unique(),columns=['Gene'])

	df_cd48_top_genes = pd.merge(df_cd48_top_genes,df_cd4[['geneSymbol','cluster.name']],how='left',right_on='geneSymbol',left_on='Gene')
	df_cd48_top_genes = pd.merge(df_cd48_top_genes,df_cd8[['geneSymbol','cluster.name']],how='left',right_on='geneSymbol',left_on='Gene')

	df_cd48_top_genes = df_cd48_top_genes[['Gene', 'cluster.name_x', 'cluster.name_y']]

	df_cd48_top_genes.columns=['Gene','CD4','CD8']
	df_cd48_top_genes['CD4'] = df_cd48_top_genes['CD4'].fillna('NotMarker')
	df_cd48_top_genes['CD8'] = df_cd48_top_genes['CD8'].fillna('NotMarker')

	df_cd48_top_genes.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'top_'+str(top_n)+'_genes_topic_CDMarkers.tsv',sep='\t',index=False)

	df_top_genes.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'top_'+str(top_n)+'_genes_topic.tsv',sep='\t',index=False)


def tcell_sample_cells_with_latent(args,cell_n=100):

	args_home = os.environ['args_home']

	df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_hh_data.tsv',sep='\t',compression='gzip')


	df_h['cluster'] = prep.cellid_to_meta_single(df_h['cell'].values,args)

	markers=['CD4.c04.Tn.il7r', \
				'CD4.c07.Tm.ANXA2', \
				'CD4.c12.Tem.GZMK', \
				'CD4.c13.Temra.CX3CR1', \
				'CD4.c14.Th17.SLC4A10', \
				'CD4.c16.Tfh.CXCR5', \
				'CD4.c17.TfhTh1.CXCL13',\
				'CD4.c20.Treg.TNFRSF9',\
				'CD4.c21.Treg.OAS1',\
				'CD4.c22.ISG.IFIT1',\
				'CD8.c01.Tn.MAL',\
				'CD8.c05.Tem.CXCR5',\
				'CD8.c06.Tem.GZMK',\
				'CD8.c07.Temra.CX3CR1',\
				'CD8.c08.Tk.TYROBP',\
				'CD8.c11.Tex.PDCD1',\
				'CD8.c15.ISG.IFIT1',\
				'CD8.c16.MAIT.SLC4A10',\
				'CD8.c17.Tm.NME1'
			   ]

	df_h = df_h[df_h['cluster'].isin(markers)]

	df_h_sample = df_h.groupby('cluster').sample(cell_n, random_state=1)

	df_h_sample.columns = [ x.replace('hh','k') for x in df_h_sample.columns]
	df_h_sample.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'hh_cell_topic_sample_CD4CD8.tsv',sep='\t',index=False)

def top_genes_marker_correlation(args,top_n=5):

	args_home = os.environ['args_home']

	df_genes=pd.read_pickle(args_home+args.input+args.raw_data_genes)

	immuneindex_file = args_home+args.input+args.nbr_model['sparse_immune_label']
	immunedat = np.load(immuneindex_file,allow_pickle=True)
	immuneindx = np.array([x[0] for x in immunedat['idx']]).astype(np.int32)
	otherindx = np.array([x for x in range(df_genes.shape[0]) if x not in immuneindx]).astype(np.int32)

	df_beta1 = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_beta1_data.tsv',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_beta2_data.tsv',sep='\t',compression='gzip')

	df_beta1.columns = df_genes.iloc[immuneindx][0].values
	df_beta2.columns = df_genes.iloc[otherindx][0].values

	db = args_home+args.database+args.tcell_signature_genes

	tag='CD8'
	df_cd4 = pd.read_excel(db,sheet_name=tag,header=1,usecols='A,B,C' )
	df_cd4 = df_cd4.sort_values('comb.ES',ascending=False)

	df_cd4_top5 = df_cd4.groupby('cluster.name').head(5).reset_index(drop=True)
	df_cd4_top5 = df_cd4_top5.sort_values('cluster.name')

	topic_clust_query = []
	clusters = df_cd4_top5['cluster.name'].unique()
	topic_clust = []
	for topic in range(df_beta1.shape[0]):
		topic_genes = df_beta1.T.iloc[:,topic].sort_values(0,ascending=False)[:top_n].index

		topic_genes_vals = df_beta1.T.iloc[:,topic].sort_values(0,ascending=False)[:top_n].values
		cell_clust = []
		for clust in clusters :
			clust_genes = df_cd4_top5[df_cd4_top5['cluster.name']==clust]['geneSymbol'].values
			clust_genes_vals = df_beta1.loc[topic,clust_genes].values
			cor = np.corrcoef(topic_genes_vals,clust_genes_vals)[0,1]
			cell_clust.append(cor)
			topic_clust_query.append(
				[topic,
				topic_genes,
				topic_genes_vals,
				clust,
				clust_genes,
				clust_genes_vals,
				cor]
			)
		topic_clust.append(cell_clust)

	df_topic_cluster_corr = pd.DataFrame(topic_clust)
	df_topic_cluster_corr.columns = clusters

	df_topic_cluster_corr.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'top_'+str(top_n)+'_topic_cluster_corr'+tag+'.tsv',sep='\t',index=False)

	pd.DataFrame(topic_clust_query).to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'top_'+str(top_n)+'_topic_cluster_corr_query'+tag+'.tsv',sep='\t',index=False)

def get_lr_pair_topic_cd48_top_genes(args,top_n=1):

	args_home = os.environ['args_home']
	l_fname = args_home+args.input+args.raw_l_data_genes
	r_fname = args_home+args.input+args.raw_r_data_genes

	ligands = pd.read_pickle(l_fname)[0].values
	receptors = pd.read_pickle(r_fname)[0].values

	df_beta1 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta1_data.tsv',sep='\t',compression='gzip')
	df_beta2 = pd.read_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta2_data.tsv',sep='\t',compression='gzip')

	df_beta1.columns = receptors
	df_beta2.columns = ligands

	top_genes = []
	top_genes = prep.generate_gene_vals(df_beta1,top_n,top_genes,'receptors')
	top_genes = prep.generate_gene_vals(df_beta2,top_n,top_genes,'ligands')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])
	top_genes = df_top_genes['Gene'].unique()

	db = args_home+args.database+args.tcell_signature_genes

	df_cd4 = pd.read_excel(db,sheet_name='CD4',header=1,usecols='A,B,C' )
	df_cd4 = df_cd4.sort_values('comb.ES',ascending=False).drop_duplicates('geneSymbol')

	df_cd8 = pd.read_excel(db,sheet_name='CD8',header=1,usecols='A,B,C' )
	df_cd8 = df_cd8.sort_values('comb.ES',ascending=False).drop_duplicates('geneSymbol')

	cd48_top_genes = list(df_cd4.geneSymbol.unique())+list(df_cd8.geneSymbol.unique())

	topic_cd48_top_genes = [ x for x in cd48_top_genes if x in top_genes]

	df_db = pd.read_csv( args_home + args.database+args.lr_db,sep='\t', usecols=['lr_pair'])

	lr_topic = []
	lr_pair = []
	for lr_p in df_db['lr_pair']:
		l,r = lr_p.split('_')[0],lr_p.split('_')[1]
		if l in df_beta2.columns and r in df_beta1.columns and \
			l in topic_cd48_top_genes and r in topic_cd48_top_genes:
			lr_topic.append((df_beta2[l]+df_beta1[r])/2)
			lr_pair.append(lr_p)
	df_lr_topic = pd.DataFrame(lr_topic)
	df_lr_topic.index=lr_pair
	df_lr_topic.to_csv(args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+'top_'+str(top_n)+'_lrpair_topic_cd4cd8_top_genes.tsv',sep='\t')

def combine_topics_tcells(args):
	
	args_home = os.environ['args_home']
	model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']

	df_h_celltype = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_hh_data.tsv',sep='\t',compression='gzip')

	df_h_state = pd.read_csv(model_file+'lrnet_interaction_states.tsv',sep='\t',compression='gzip')

	df_h_state['state'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_h_state.iterrows()]

	df_h_celltype['topic'] = df_h_celltype.iloc[:,1:].idxmax(axis=1)

	dflatent = pd.merge(df_h_state[['cell','state']],df_h_celltype[['cell','topic']],how='left',on='cell')

	meta_path = args_home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep='\t')
	dflatent['label'] = pd.merge(dflatent,df_meta,right_on='cellID',left_on='cell',how='left')['meta.cluster'].values

	dfsummary = dflatent.groupby(['label','topic','state']).count()
	dfsummary = dfsummary.reset_index()
	dfsummary.to_csv(model_file+'model_summary.csv',index=False)

def cell_interaction_network(args):
	
	import igraph
	import seaborn as sns

	args_home = os.environ['args_home']

	model_file = args_home+args.output+args.lr_model['out']+'10_ld/'+args.lr_model['mfile']

	df = pd.read_csv(model_file+'lrnet_interaction_states.tsv',sep='\t',compression='gzip')

	meta_path = args_home+ args.database +args.metadata
	df_meta = pd.read_csv(meta_path,sep='\t')
	df['label'] = pd.merge(df,df_meta,right_on='cellID',left_on='cell',how='left')['meta.cluster'].values

	df['state'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df.iterrows()]

	df = df.dropna()
	df = df[df.label.str.contains('CD4')]
	df = df[df.label.str.contains('Tn|Tm|Temra|Th|Treg')]

	df = df[df.state != 9]
	df = df[df.state != 6]

	dfs = df.groupby('label', group_keys=False).apply(pd.DataFrame.sample,frac=0.1)
	cells = [x[5:10] for x in list(dfs.label.values)]

	states = ['t'+str(x) for x in dfs['state'].unique()]

	## construct graph
	if len(states)<=10:
		colors = sns.color_palette("Paired")
	elif len(states)>10:
		colors = sns.color_palette('Spectral',dfs['state'].unique().max()+1)

	g = igraph.Graph()
	g.add_vertices(cells+states)

	gedges = []
	for i,row in dfs.iterrows():
		cell = row.label[5:10]
		s = row['state']
		interaction_state = 't'+str(s)
		epair = cell + '_' + interaction_state

		if epair not in gedges:
			g.add_edge(cell,interaction_state,weight=0.1,color=colors[s])
			gedges.append(epair)
		else:
			eid = g.get_eid(cell,interaction_state)
			g.es[eid]['weight'] += 0.1


	to_delete_eids = [e.index for e in g.es if e['weight']<0.5 ]
	g.delete_edges(to_delete_eids)

	to_delete_vids = [v.index for v in g.vs if v.degree()<1]
	g.delete_vertices(to_delete_vids)

	
	for v in g.vs:
		if "." in v['name']: 
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

