
##################################################################
    # setup
##################################################################
import experiment_interaction_topic 
import pandas as pd
from pathlib import Path
from analysis import _topics,_network

server = Path.home().as_posix()
experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'
experiment_home = server+experiment_home
spr = experiment_interaction_topic.get_experiment_model(experiment_home)
print(spr.cell_topic.model_id)
print(spr.interaction_topic.model_id)

def argmax_latent_plot():
    ##################################################################
    # 1 analysis of latent h
    ##################################################################
    from analysis import _topics
    import matplotlib.pylab as plt
    plt.rcParams['figure.figsize'] = [7, 5]
    plt.rcParams['figure.autolayout'] = True
    import colorcet as cc
    import seaborn as sns

    '''
    The latent h of interaction topic is based on cell-pairs. So, we have to take argmax for each pair. 
    Each row is one cell and argmax h for its 159 neighbours.
    '''
    df_its = spr.interaction_topic_states()

    ###save/load
    df_its.to_csv(spr.interaction_topic.model_id+'_it_1_h_argmax.csv.gz',index=False,compression='gzip')
    df_its = pd.read_csv(spr.interaction_topic.model_id+'_it_1_h_argmax.csv.gz',compression='gzip')

    df_hmax = pd.DataFrame(pd.Series(df_its.iloc[:,1:].values.flatten()).value_counts()).reset_index().rename(columns={'index':'cell_topic',0:'argmax_count'})
    df_hmax = df_hmax.sort_values('argmax_count',ascending=False)
    p = sns.barplot(x='cell_topic',y='argmax_count',data=df_hmax,color='blue')
    p.set_xlabel("Topic",fontsize=20)
    p.set_ylabel("Count(argmax)",fontsize=20)
    plt.savefig(spr.interaction_topic.model_id+'_it_1_h_argmax.png');plt.close()

def metafile():
    ##################################################################
    # make metafile
    ##################################################################

    '''
    get cell ids and cell topic from cell topic model
    and get argmax cell pair neighbour toopic 
    and max proportion interaction topic for each cell pair
    '''
    df_kmeans = pd.read_csv(spr.cell_topic.model_id +'_ct_3_h_kmeans.csv.gz')
    df_its = pd.read_csv(spr.interaction_topic.model_id+'_it_1_h_argmax.csv.gz',compression='gzip')

    df_its['interact_topic'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_its.iterrows()]
    df_kmeans = pd.merge(df_kmeans,df_its[['cell','interact_topic']],how='left',on='cell')
    # df_kmeans = df_kmeans.drop(columns=['umap1', 'umap2', 'label', 'l1','celltype_minor','celltype',])
    df_kmeans.to_csv(spr.interaction_topic.model_id+'_it_2_meta_ct_argmax_maxprop.csv.gz',index=False,compression='gzip')

def metafile_nbrs():
    df_kmeans = pd.read_csv(spr.cell_topic.model_id +'_ct_3_h_kmeans.csv.gz')
    df_its = pd.read_csv(spr.interaction_topic.model_id+'_it_1_h_argmax.csv.gz',compression='gzip')
    df_nbr = spr.cell_topic.neighbour.copy()

    selected_ctypes=['Cancer','Epithelial','T','B_n|B_Mem','Monocyte','Endothelial','PVL','CAF']
    for selected_ct in selected_ctypes:
        selected_cells = df_kmeans[df_kmeans['cluster_celltype'].str.contains(selected_ct)]['cell'].values
        res = _topics.get_cell_neighbours_states(df_nbr,df_its,selected_cells)
        df_res = pd.DataFrame(res,columns=[selected_ct,'nbr','interact_topic'])
        df_res = df_res.explode(['nbr','interact_topic'])
        df_res.to_csv(spr.interaction_topic.model_id+'_it_2_meta_'+selected_ct+'_cells_nbrs.csv.gz',index=False,compression='gzip')

def celltype_ct_it_distribution():

    ##################################################################
    # cell types from cell topic and 
    # sub types from interaction topic analysis
    # ** cancer cells summary plot **
    ##################################################################

    selected_ctypes=['Cancer','Epithelial','T','B_n|B_Mem','Monocyte','Endothelial','PVL','CAF']

    df_combined = pd.DataFrame()
    for celltype in selected_ctypes:
        print(celltype)
        df=pd.read_csv(spr.interaction_topic.model_id+'_it_2_meta_'+celltype+'_cells_nbrs.csv.gz',compression='gzip')
        selected_int_topics = [2,4,7,10,18,22,24]
        df = df[df['interact_topic'].isin(selected_int_topics)]
        dfmeta = pd.read_csv(spr.interaction_topic.model_id+'_it_2_meta_ct_argmax_maxprop.csv.gz')
        df = pd.merge(df,dfmeta[['cell','cell_topic']],left_on=celltype,right_on='cell',how='left')
        dfs = df.groupby(['cell_topic','interact_topic'])[celltype].count().reset_index()
        df_select_topic  = pd.DataFrame(dfmeta[dfmeta['cluster_celltype'].str.contains(celltype)]['cell_topic'].value_counts())
        selected_cell_topic=list(df_select_topic[df_select_topic['cell_topic']>100].index)
        dfs = dfs[dfs['cell_topic'].isin(selected_cell_topic)]
        celltopic_sum = dict(dfs.groupby('cell_topic')[celltype].sum())
        dfs['ncount'] = [x/celltopic_sum[y] for x,y in zip(dfs[celltype],dfs['cell_topic'])]

        dfs = dfs.rename(columns={celltype:'celltype_val'})
        dfs['celltype'] = celltype

        df_combined = pd.concat([df_combined, dfs], axis=0, ignore_index=True)

    df_combined.to_csv(spr.interaction_topic.model_id+'_it_3_cell_nbr_ct_it_summary.csv.gz',index=False,compression='gzip')

def celltype_distribution_interaction_topic():
    ##################################################################
    # 3.3 
    # cell type distribution of interaction topic of cancer cells
    # it with interesting pattern - 2,4,7,10,18,22,24

    df=pd.read_csv(spr.interaction_topic.model_id+'_it_2_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')

    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]
    dfmeta = pd.read_csv(spr.interaction_topic.model_id+'_it_2_meta_ct_argmax_maxprop.csv.gz')
    df = pd.merge(df,dfmeta[['cell','cluster_celltype']],left_on='nbr',right_on='cell',how='left')

    df = df[~df['cluster_celltype'].str.contains('Cancer')]

    df_grp = df.groupby(['cluster_celltype','interact_topic'])['interact_topic'].size().rename('count').reset_index()
    ##normalize
    celltype_sum = dict(df_grp.groupby('interact_topic')['count'].sum())

    df_grp['ncount'] = [x/celltype_sum[y] for x,y in zip(df_grp['count'],df_grp['interact_topic'])]

    df_grp.to_csv(spr.interaction_topic.model_id+'_it_4_cancercells_it_celltypedist.csv.gz',index=False,compression='gzip')

def cluster_celltype_distribution_interaction_topic():
    ##################################################################
    # 3.3 
    # cell type distribution of interaction topic of cancer cells
    # check for normal datasets bias

    df=pd.read_csv(spr.interaction_topic.model_id+'_it_2_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')

    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]
    dfmeta = pd.read_csv(spr.interaction_topic.model_id+'_it_2_meta_ct_argmax_maxprop.csv.gz')

    dfmeta['cluster_celltype'] = [str(x)+'/'+str(y) for x,y in zip(dfmeta['cluster'],dfmeta['cluster_celltype'])]

    df = pd.merge(df,dfmeta[['cell','cluster_celltype']],left_on='nbr',right_on='cell',how='left')

    df = df[~df['cluster_celltype'].str.contains('Cancer')]

    df_grp = df.groupby(['cluster_celltype','interact_topic'])['interact_topic'].size().rename('count').reset_index()
    ##normalize

    ### check for normal and not normal
    df_grp = df_grp[~df_grp['cluster_celltype'].str.contains('_n')]

    celltype_sum = dict(df_grp.groupby('interact_topic')['count'].sum())

    df_grp['ncount'] = [x/celltype_sum[y] for x,y in zip(df_grp['count'],df_grp['interact_topic'])]

    
    df_grp.to_csv(spr.interaction_topic.model_id+'_it_4_cancercells_it_cluster_celltypedist_notnormal.csv.gz',index=False,compression='gzip')

def raw_topgenes_exp_heatmap():
    # ##################################################################
    # # 3.2 
    # # get mean expression of cancer nbr cells from above sample dataset
    # # draw heatmap
    # '''
    # -for each interaction topic get cancer cells select its neighbour cells
    # -get raw data for these cancer cells and nbr cells
    # - normalize rowwise such that each ligand expression is proportion of all ligand expression
    # - get top 25 ligand and receptor genes for interaction topic
    # - filter those top genes and take mean, add 1e-5, and apply sqrt for heatmap
    # '''

    df=pd.read_csv(spr.interaction_topic.model_id+'_it_2_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')
    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]
    
    _topics.plt_cn(spr,df,selected_int_topics)

def top_genes():
    ##################################################################
    # 2 top genes
    ##################################################################


    from analysis import _topics

    #### get top genes based on beta weight matrix - l/r separately

    top_n = 10
    df_top_genes = _topics.topic_top_lr_genes(spr,top_n)
    df_top_genes.to_csv(spr.interaction_topic.model_id+'_it_5_beta_weight_top_'+str(top_n)+'_genes.csv.gz',index=False,compression='gzip')

def correlation_lr_network():
    ##################################################################
    # 4 gene gene L/R correlation network 
    ##################################################################

    ##################################################################
    '''
    For each interaction topics,
    for cells -> combine cells from all cell pairs 
    for genes -> combine ligand receptor genes
    construct gene correlation matrix using raw count data as gene expression, add 1 where zero.
    '''
    from analysis import _network
    from util._io import read_config
    from collections import namedtuple
    from analysis import _topics


    # first get neighbour cells of each cancer cell and their interaction topic 
    '''
    total data size will be 4185039
    total cancer cells 26321
            26321 x 159 neighbours = 4185039
    '''

    df=pd.read_csv(spr.interaction_topic.model_id+'_it_2_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')

    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]
    dfmeta = pd.read_csv(spr.interaction_topic.model_id+'_it_2_meta_ct_argmax_maxprop.csv.gz')
    df = pd.merge(df,dfmeta[['cell','cluster_celltype']],left_on='nbr',right_on='cell',how='left')

    experiment_config = read_config(experiment_home+'config.yaml')
    args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())
    spr.data.raw_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_data)

    corr_th = 0.25
    zcutoff = 0.5
    selected_int_topics = [2,4,7,10,18,22,24]
    _network.lr_correlation_network(spr,df,selected_int_topics,zcutoff,corr_th,'lr')
    _network.lr_correlation_network(spr,df,selected_int_topics,zcutoff,corr_th,'ll')
    _network.lr_correlation_network(spr,df,selected_int_topics,zcutoff,corr_th,'rr')

def struct_plot():
    ##################################################################
    # 3 interaction topic struct plot 
    ##################################################################

    ##################################################################
    # 3.1 
    # get samples of cell pairs interaction topic
    '''
    take 100 cells sample from each cancer cell topic with >100 cells
    in total 13 cell topics
    get 159 neighbours and interaction topic probability for each cancer cell
    combine data 
    100 cancer cells sample from 13 cell topics is 1300

    1300 * 159 = 206700 cell pairs

    using kmeans identify interaction topic dominated cluster id and take 50 samples from each cluster for structure plot
    '''

    import numpy as np
    from sklearn.cluster import KMeans

    df_kmeans = pd.read_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_kmeans.csv.gz')
    df_kmeans = df_kmeans[ df_kmeans.cluster_celltype.str.contains('Cancer')]

    # identify cancer cell topic with > 100 cells
    df_kmeans = df_kmeans[ df_kmeans['cell_topic'].isin(df_kmeans.topic.value_counts().index[:13])]
    df_kmeans = df_kmeans.groupby('cell_topic').sample(n=100, random_state=1)
    cancer_cells = df_kmeans['cell'].values

    dfh_sample = spr.interaction_topic_prop_with_cellids(cancer_cells)

    kmeans = KMeans(n_clusters=100, random_state=0).fit(dfh_sample.iloc[:,:-2].to_numpy())
    dfh_sample['cluster'] = kmeans.labels_

    ###save
    dfh_sample.to_csv(spr.interaction_topic.model_id+'_it_topic_sample_kmeans.csv.gz',index=False,compression='gzip')


    ### plot entire kmeans cluster
    dfh_sample = pd.read_csv(spr.interaction_topic.model_id+'_it_topic_sample_kmeans.csv.gz')
    dfh_sample_kmeans = dfh_sample.groupby('cluster').sample(n=50, random_state=1)
    dfh_sample_kmeans['cell'] = [ x+'/'+y for x,y in  zip(dfh_sample_kmeans['cancer_cell'],dfh_sample_kmeans['nbr'])]
    dfh_sample_kmeans = dfh_sample_kmeans.drop(columns=['cancer_cell','nbr'])
    dfh_sample_kmeans = dfh_sample_kmeans.rename(columns={'cluster':'Topic'})
    dfh_sample_kmeans.to_csv(spr.interaction_topic.model_id+'_it_h_sample_kmeans.csv.gz',index=False,compression='gzip')

    ### plot interaction topic enriched kmeans cluster
    '''
    enriched topics and cluster pairs are identified manually from the above 
    entire kmeans cluster plot
    '''
    dfh_sample = pd.read_csv(spr.interaction_topic.model_id+'_it_topic_sample_kmeans.csv.gz')
    dfh_sample_kmeans = dfh_sample.groupby('cluster').sample(n=50, random_state=1)
    dfh_sample_kmeans['cell'] = [ x+'/'+y for x,y in  zip(dfh_sample_kmeans['cancer_cell'],dfh_sample_kmeans['nbr'])]
    dfh_sample_kmeans = dfh_sample_kmeans.drop(columns=['cancer_cell','nbr'])
    dfh_sample_kmeans = dfh_sample_kmeans.rename(columns={'cluster':'Topic'})

    clust_it_pair = [33,19,4,97,72,67,5,75,52,79,90,98,6,43,93,83,55,0,68,82,86,12,10,69,57 ]
    dfh_sample_kmeans = dfh_sample_kmeans[dfh_sample_kmeans['Topic'].isin(clust_it_pair)]
    dfh_sample_kmeans['Topic'] = ['t'+str(i) for i in dfh_sample_kmeans['Topic']]

    torder={
    't33':1,
    't19':2,
    't4':3,
    't97':4,
    't72':5,
    't67':6,
    't5':7,
    't75':8,
    't52':9,
    't79':10,
    't90':11,
    't98':12,
    't6':13,
    't43':14,
    't93':15,
    't83':16,
    't55':17,
    't0':18,
    't68':19,
    't82':20,
    't86':21,
    't12':22,
    't10':23,
    't69':24,
    't57':25
    }

    dfh_sample_kmeans['Topic'] = [torder[i] for i in dfh_sample_kmeans['Topic']]

    selected_int_topics = [3,11,13,16,18,20,23]

    dfh_sample_kmeans = dfh_sample_kmeans[dfh_sample_kmeans['Topic'].isin(selected_int_topics)]

    dfh_sample_kmeans = dfh_sample_kmeans.sort_values('Topic')
    dfh_sample_kmeans.to_csv(spr.interaction_topic.model_id+'_it_h_sample_kmeans_selected.csv.gz',index=False,compression='gzip')

def deg_analysis():
    ##################################################################
    # DEG analysis
    ##################################################################

    selected_int_topics = [2,4,10,18,22,24]

    df = pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics.csv.gz',compression='gzip')
    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]

    df_its = pd.read_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz',compression='gzip')
    df_nbr = spr.cell_topic.neighbour

    res = _topics.get_cell_neighbours_states(df_nbr,df_its,df['cell'].values)
    df_res = pd.DataFrame(res,columns=['cell_id','nbr','interact_topic'])
    df_res = df_res.explode(['nbr','interact_topic'])
    df_res = pd.merge(df_res,df[['cell','cell_topic','cluster','cluster_celltype']],left_on='nbr',right_on='cell',how='left')
    df_res = df_res.drop(columns=['cell'])

    df_res = df_res[df_res['interact_topic'].isin(selected_int_topics)]

    df_res.to_pickle(spr.interaction_topic.model_id+'_it_model_all_topics_all_cells.pkl')

    df_res = pd.read_pickle(spr.interaction_topic.model_id+'_it_model_all_topics_all_cells.pkl')

    query_cells = df_res['cell_id'].values
    query_cells = pd.Series(query_cells).unique()
    df = spr.interaction_topic_prop_with_cellids_nbrsummed(query_cells)
    df.to_csv(spr.interaction_topic.model_id+'_it_model_deg_analysis.csv.gz',index=False,compression='gzip')

    #### DEG data check
    df = pd.read_csv(spr.interaction_topic.model_id+'_it_model_deg_analysis.csv.gz')
    df_meta = pd.read_csv(spr.interaction_topic.model_id+'_alltopics_meta.csv.gz',compression='gzip')
    dftest = df_meta[df_meta['cell'].isin(df['cell_id'].values)]
    dftest.shape
    dftest['interact_topic'].value_counts()
    dftest.cluster_celltype.value_counts()

# ##################################################################
# ##################################################################
# # NOT INCLUDED IN FIGURES
# ##################################################################

# ##################################################################
# # OTHER cancer centric view
# ##################################################################

# ##################################################################
# # 4.1 
# # generate topics summary 

# '''
# Take umap coordinates, topic assignment, cluster assignment, cell type, and cluster majority cell type from
# cell topic analysis and merge with state assignment from interaction topic.
# '''
# spr.interaction_topic.neighbour_h =  pd.read_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz')
# df_kmeans = pd.read_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_kmeans.csv.gz')

# df_all_topics = _topics.get_topics(spr,df_kmeans)
# df_all_topics.to_csv(spr.interaction_topic.model_id+'_it_model_all_topics.csv.gz',index=False,compression='gzip')


# ##################################################################
# # 5 network graph
# # Check interaction state of cancer cells and neighbouring cells 

# '''
# take all cancer cells and remove state 7 and 10 
# plot network to show cancer cell states and neighbour cells belonging
# to that state
# '''
# from analysis import _network
# df = pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics_cancer_cells.csv.gz',compression='gzip')
# # df = df[df['cluster_celltype'] == 'Cancer Epithelial']
# # df = df[~df['interact_topic'].isin([7,10])]
# _network.cell_interaction_network(spr,df)


# ##################################################################
# # 6 topic wise lr network graph
# ##################################################################

# from analysis import _network

# df_db = pd.read_csv( experiment_home + spr.args.database+ spr.args.lr_db,sep='\t', usecols=['lr_pair'])

# states = [10,7,4,2,22,24,18,1]
# top_lr = 200
# keep_db= True
# _network.interaction_statewise_lr_network(spr,states,top_lr,keep_db,df_db)

# top_lr = 25
# keep_db=False
# _network.interaction_statewise_lr_network(spr,states,top_lr,keep_db,df_db)


# ##################################################################
# # 7 interaction topic gse analysis
# ##################################################################
# '''
# take ligand and receptor weight and average them
# '''
# df_db = pd.read_csv( experiment_home + spr.args.database+ spr.args.lr_db,sep='\t', usecols=['lr_pair'])
# df = _gsea.gse_interactiontopic_v2(spr,df_db)
# df.to_csv(spr.interaction_topic.model_id+'_it_gsea.csv.gz',index=False,compression='gzip')

# from analysis import _gsea
# _gsea.gsea(spr)


