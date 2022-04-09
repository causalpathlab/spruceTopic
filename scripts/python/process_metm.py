import os
import pandas as pd
import numpy as np
from collections import namedtuple
from gen_util.io import read_config

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
    
    immuneindex_file = args_home+args.input+args.nbr_model['sparse_immune_label']
    immunedat = np.load(immuneindex_file,allow_pickle=True)
    immuneindx = np.array([x[0] for x in immunedat['idx']]).astype(np.int32)
    otherindx = np.array([x for x in range(df_genes.shape[0]) if x not in immuneindx]).astype(np.int32)

    df_beta1 = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_beta1_data.tsv',sep='\t',compression='gzip')
    df_beta2 = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_beta2_data.tsv',sep='\t',compression='gzip')

    df_beta1.columns = df_genes.iloc[immuneindx][0].values
    df_beta2.columns = df_genes.iloc[otherindx][0].values

    top_genes = []
    top_genes = generate_gene_vals(df_beta1,top_n,top_genes,'immune')
    top_genes = generate_gene_vals(df_beta2,top_n,top_genes,'non-immune')

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

def sample_cells_with_latent(args,cell_n=50):

    from sklearn.cluster import KMeans

    args_home = os.environ['args_home'] 

    df_h = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'etm_hh_data.tsv',sep='\t',compression='gzip')
    kmeans = KMeans(n_clusters=df_h.shape[1]-1, random_state=0).fit(df_h.iloc[:,1:].to_numpy())
    df_h['cluster'] = kmeans.labels_
    df_h_sample = df_h.groupby('cluster').sample(cell_n, random_state=1)

    # df_h_sample = pd.DataFrame
    # for i in range(16):
    #     df_h = df_h.sort_values('hh'+str(i))
    #     if i ==0:
    #         df_h_sample = df_h.iloc[:10,:]
    #     else:
    #         df_h_sample = pd.concat([df_h_sample,df_h.iloc[:10,:]],axis=0, ignore_index=True)
    # df_h_sample = df_h_sample.drop_duplicates()
    df_h_sample.columns = [ x.replace('hh','') for x in df_h_sample.columns]
    df_h_sample.to_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'hh_cell_topic_sample.tsv',sep='\t',index=False)

def get_top_lr_genes(args,top_n=3):


	df_genes=pd.read_pickle(args.home+args.input+args.raw_genes)
	
	fname = args.home+args.input+args.raw_data
	l_fname = fname.replace('.pkl','.ligands.pkl')
	r_fname = fname.replace('.pkl','.receptors.pkl')


	df_l = pd.read_pickle(l_fname)
	df_r = pd.read_pickle(r_fname)
	ligands = df_l.columns[1:]
	receptors = df_r.columns[1:]

	df_beta1 = pd.read_csv(args.home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta1_data.csv',sep='\t')
	df_beta2 = pd.read_csv(args.home+args.output+args.lr_model['out']+args.lr_model['mfile']+'etm_beta2_data.csv',sep='\t')

	df_beta1.columns = receptors
	df_beta2.columns = ligands

	top_genes = []
	top_genes = generate_gene_vals(df_beta1,top_n,top_genes,'receptors')
	top_genes = generate_gene_vals(df_beta2,top_n,top_genes,'ligands')

	df_top_genes = pd.DataFrame(top_genes,columns=['Topic','GeneType','Genes','Gene','Proportion'])
	df_top_genes.to_csv(args.home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'top_'+str(top_n)+'_genes_topic.tsv',sep='\t',index=False)

    