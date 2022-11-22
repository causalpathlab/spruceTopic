
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
print(spr.cell_topic.id)
print(spr.interaction_topic.id)

def argmax_latent_plot():
    ##################################################################
    # 1 analysis of latent h
    ##################################################################
    import matplotlib.pylab as plt
    plt.rcParams['figure.figsize'] = [20, 5]
    plt.rcParams['figure.autolayout'] = True
    import seaborn as sns

    '''
    The latent h of interaction topic is based on cell-pairs. So, we have to take argmax for each pair. 
    Each row is one cell and argmax h for its 159 neighbours.
    '''
    df_its = spr.interaction_topic_states()
    df_its.to_csv(spr.interaction_topic.
    id+'1_h_argmax.csv.gz',index=False,compression='gzip')

    df_its = pd.read_csv(spr.interaction_topic.id+'1_h_argmax.csv.gz',compression='gzip')

    df_hmax = pd.DataFrame(pd.Series(df_its.iloc[:,1:].values.flatten()).value_counts()).reset_index().rename(columns={'index':'cell_topic',0:'argmax_count'})
    df_hmax = df_hmax.sort_values('argmax_count',ascending=False)
    p = sns.barplot(x='cell_topic',y='argmax_count',data=df_hmax,palette=sns.color_palette("hls", 25))
    p.set_xlabel("Interaction topic",fontsize=20)
    p.set_ylabel("Total cell pair",fontsize=20)
    plt.savefig(spr.interaction_topic.id+'1_h_argmax.pdf',dpi=300);plt.close()

def top_genes_tpwise():
    top_n = 10
    df_topl = _topics.generate_top_genes_topicwise(spr.interaction_topic.beta_lm,top_n)
    df_topl['type']='ligand'
    df_topr = _topics.generate_top_genes_topicwise(spr.interaction_topic.beta_rm,top_n)
    df_topr['type']='receptor'
    df_top = pd.concat([df_topl,df_topr], axis=0, ignore_index=True)
    df_top.to_csv(spr.interaction_topic.id+'2_beta_weight_top_'+str(top_n)+'_genes_topicwise.csv.gz',index=False,compression='gzip')

    dfl,dfr = _topics.get_zscores(spr)
    dfl.to_csv(spr.interaction_topic.id+'2_zscore_ligands.csv',index=False)
    dfr.to_csv(spr.interaction_topic.id+'2_zscore_receptors.csv',index=False)

def top_genes():
    ##################################################################
    # 2 top genes
    ##################################################################


    from analysis import _topics

    #### get top genes based on beta weight matrix - l/r separately

    top_n = 10
    df_top_genes = _topics.topic_top_lr_genes(spr,top_n)
    df_top_genes.to_csv(spr.interaction_topic.id+'2_beta_weight_top_'+str(top_n)+'_genes.csv.gz',index=False,compression='gzip')
    
    top_n = 600
    df_top_genes = _topics.topic_top_lr_genes(spr,top_n)
    df_top_genes.to_csv(spr.interaction_topic.id+'2_beta_weight_top_'+str(top_n)+'_genes.csv.gz',index=False,compression='gzip')

def metafile():
    ##################################################################
    # make metafile
    ##################################################################

    '''
    get cell ids and cell topic from cell topic model
    and get argmax cell pair neighbour toopic 
    and max proportion interaction topic for each cell pair
    '''
    df_kmeans = pd.read_csv(spr.cell_topic.id +'2_i_kmeans.csv.gz')
    df_its = pd.read_csv(spr.interaction_topic.id+'1_h_argmax.csv.gz',compression='gzip')

    df_its['interact_topic'] = [ pd.Series(vals).value_counts().index[0] for indx,vals in df_its.iterrows()]
    df_kmeans = pd.merge(df_kmeans,df_its[['cell','interact_topic']],how='left',on='cell')
    # df_kmeans = df_kmeans.drop(columns=['umap1', 'umap2', 'label', 'l1','celltype_minor','celltype',])
    xrep={
        
        'B_Memory':'B', 
        'B_n':'B', 
        
        'CAFs_MSC_iCAF_like':'CAF', 
        'CAFs_myCAF_like':'CAF',
        'CAFs_n':'CAF', 
        
        'Cancer_Basal_SC':'Cancer', 
        'Cancer_Cycling':'Cancer',
        'Cancer_Her2_SC':'Cancer',
        'Cancer_LumA_SC':'Cancer',
        'Cancer_LumB_SC':'Cancer',

        'Endothelial_ACKR1':'Endothelial',
        'Endothelial_n':'Endothelial', 
        
        'Epithelial_Luminal':'Epithelial',
        'Epithelial_Myoepi':'Epithelial',
        'Epithelial_n':'Epithelial',
        
        'Macrophage':'Myeloid', 
        'Macrophage_n':'Myeloid',  
        'Monocyte':'Myeloid', 

        'PVL_Differentiated':'PVL',
        'PVL_Immature':'PVL',
        
        'Plasma':'Plasma',
        'Plasma_n':'Plasma',
        
        'T_CD4':'T', 
        'T_CD8_pan':'T'

    }

    df_kmeans['cluster_celltype'] = [ xrep[x] for x in df_kmeans['cluster_celltype']]

    df_kmeans.to_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz',index=False,compression='gzip')

def metafile_nbrs():
    df_kmeans = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')

    df_its = pd.read_csv(spr.interaction_topic.id+'1_h_argmax.csv.gz',compression='gzip')
    df_nbr = spr.cell_topic.neighbour.copy()

    selected_ctypes=['T', 'Cancer', 'B', 'Epithelial', 'Myeloid', 'Endothelial', 'Plasma',
       'CAF', 'PVL']
    for selected_ct in selected_ctypes:
        print(selected_ct)
        selected_cells = df_kmeans[df_kmeans['cluster_celltype']==selected_ct]['cell'].values
        res = _topics.get_cell_neighbours_states(df_nbr,df_its,selected_cells)
        df_res = pd.DataFrame(res,columns=[selected_ct,'nbr','interact_topic'])
        df_res = df_res.explode(['nbr','interact_topic'])
        df_res.to_csv(spr.interaction_topic.id+'3_meta_'+selected_ct+'_cells_nbrs.csv.gz',index=False,compression='gzip')

def celltype_it_distribution():

    ##################################################################
    # cell types from cell topic and 
    # sub types from interaction topic analysis
    # ** cancer cells summary plot **
    ##################################################################

    selected_ctypes=['T', 'Cancer', 'B', 'Epithelial', 'Myeloid', 'Endothelial', 'Plasma',
       'CAF', 'PVL']

    df_combined = pd.DataFrame()
    for celltype in selected_ctypes:
        print(celltype)
        df=pd.read_csv(spr.interaction_topic.id+'3_meta_'+celltype+'_cells_nbrs.csv.gz',compression='gzip')
        selected_int_topics = [2,4,7,10,18,22,24]
        df = df[df['interact_topic'].isin(selected_int_topics)]
        dfs = df.groupby(['interact_topic'])[celltype].count().reset_index()
        itopic_sum = dfs[celltype].sum()
        dfs['ncount'] = [x/itopic_sum for x in dfs[celltype]]

        dfs = dfs.rename(columns={celltype:'celltype_val'})
        dfs['celltype'] = celltype

        df_combined = pd.concat([df_combined, dfs], axis=0, ignore_index=True)

    df_combined.to_csv(spr.interaction_topic.id+'4_cell_nbr_it_summary.csv.gz',index=False,compression='gzip')

def celltype_ct_it_distribution():

    ##################################################################
    # cell types from cell topic and 
    # sub types from interaction topic analysis
    # ** cancer cells summary plot **
    ##################################################################

    selected_ctypes=['T', 'Cancer', 'B', 'Epithelial', 'Myeloid', 'Endothelial', 'Plasma',
       'CAF', 'PVL']

    df_combined = pd.DataFrame()
    for celltype in selected_ctypes:
        print(celltype)
        df=pd.read_csv(spr.interaction_topic.id+'3_meta_'+celltype+'_cells_nbrs.csv.gz',compression='gzip')
        selected_int_topics = [2,4,7,10,18,22,24]
        df = df[df['interact_topic'].isin(selected_int_topics)]
        dfmeta = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')
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

    df_combined.to_csv(spr.interaction_topic.id+'4_cell_nbr_ct_it_summary.csv.gz',index=False,compression='gzip')

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

    df=pd.read_csv(spr.interaction_topic.id+'3_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')
    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]
        
    _topics.plt_cn(spr,df,selected_int_topics)

def it_celltypedist():
    ##################################################################
    # 3.3 
    # cell type distribution of interaction topic of cancer cells
    # it with interesting pattern - 2,4,7,10,18,22,24

    selected_ctypes=['T', 'Cancer', 'B', 'Epithelial', 'Myeloid', 'Endothelial', 'Plasma',
       'CAF', 'PVL']
    
    df_combined = pd.DataFrame()

    for celltype in selected_ctypes:

        print(celltype)
        
        df=pd.read_csv(spr.interaction_topic.id+'3_meta_'+celltype+'_cells_nbrs.csv.gz',compression='gzip')

        selected_int_topics = [2,4,7,10,18,22,24]
        df = df[df['interact_topic'].isin(selected_int_topics)]
        dfmeta = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')
        df = pd.merge(df,dfmeta[['cell','cluster_celltype']],left_on='nbr',right_on='cell',how='left')

        # df = df[~df['cluster_celltype'].str.contains('Cancer')]

        df_grp = df.groupby(['cluster_celltype','interact_topic'])['interact_topic'].size().rename('count').reset_index()
        ##normalize
        celltype_sum = dict(df_grp.groupby('interact_topic')['count'].sum())

        df_grp['ncount'] = [x/celltype_sum[y] for x,y in zip(df_grp['count'],df_grp['interact_topic'])]

        df_grp = df_grp.sort_values('cluster_celltype')

        df_grp['celltype'] = celltype

        df_combined = pd.concat([df_combined, df_grp], axis=0, ignore_index=True)

    df_combined.to_csv(spr.interaction_topic.id+'4_it_celltypedist.csv.gz',index=False,compression='gzip')
    
def struct_plot():
    import numpy as np
    from sklearn.cluster import KMeans

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


    df=pd.read_csv(spr.interaction_topic.id+'3_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')

    selected_int_topics = [2,4,7,10,18,22,24]

    df = df[df['interact_topic'].isin(selected_int_topics)]

    df = df.groupby(['Cancer','interact_topic']).sample(n=1, random_state=1)

    cancer_nbr_cells = [[x,y] for x,y in zip(df['Cancer'],df['nbr'])]
    
    tp_prob = spr.interaction_topic_prop_with_cellids(cancer_nbr_cells)

    df_sample = pd.DataFrame([x[2][0] for x in tp_prob])
    df_sample['cancer_indx'] = [x[0] for x in tp_prob]
    df_sample['nbr_indx'] = [x[1] for x in tp_prob]


    df_it = df_sample[ [x for x  in selected_int_topics]]
    df_it =  df_it.div(df_it.sum(axis=1), axis=0)
    df_it[['cancer_indx','nbr_indx']] = df_sample[['cancer_indx','nbr_indx']]
    

    kmeans = KMeans(n_clusters=7, random_state=0).fit(df_it.iloc[:,:-2].to_numpy())
    df_it['cluster'] = kmeans.labels_

    df_it['cell'] = [ str(x)+'_'+str(y) for x,y in  zip(df_it['cancer_indx'],df_it['nbr_indx'])]

    df_it = df_it.drop(columns=['cancer_indx','nbr_indx'])
    df_it = df_it.rename(columns={'cluster':'Topic'})

    df_it.to_csv(spr.interaction_topic.id+'6_topic_sample_kmeans_stplot.csv.gz',index=False,compression='gzip')

def caner_celltype_it_distribution():

    ##################################################################
    # cell types from cell topic and 
    # sub types from interaction topic analysis
    # ** cancer cells summary plot **
    ##################################################################

    selected_stypes=['TNBC', 'ER+', 'HER2+']
    celltype = 'Cancer'
    df_combined = pd.DataFrame()
    for subtype in selected_stypes:
        print(subtype)
        df=pd.read_csv(spr.interaction_topic.id+'3_meta_'+celltype+'_cells_nbrs.csv.gz',compression='gzip')
        selected_int_topics = [2,4,7,10,18,22,24]
        df = df[df['interact_topic'].isin(selected_int_topics)]

        dfmeta = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')

        df = pd.merge(df,dfmeta[['cell','subtype']],left_on=celltype,right_on='cell',how='left')

        df = df[df['subtype']==subtype]

        dfs = df.groupby(['interact_topic'])[celltype].count().reset_index()
        itopic_sum = dfs[celltype].sum()
        dfs['ncount'] = [x/itopic_sum for x in dfs[celltype]]

        dfs = dfs.rename(columns={celltype:'celltype_val'})
        dfs['celltype'] = subtype

        df_combined = pd.concat([df_combined, dfs], axis=0, ignore_index=True)

    df_combined.to_csv(spr.interaction_topic.id+'5_cancer_stype_nbr_it_summary.csv.gz',index=False,compression='gzip')

def cancer_subtype_ct_it_distribution():

    selected_stypes=['TNBC', 'ER+', 'HER2+']
    df_combined = pd.DataFrame()
    celltype='Cancer'

    for subtype in selected_stypes:
        print(subtype)
        df=pd.read_csv(spr.interaction_topic.id+'3_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')
        selected_int_topics = [2,4,7,10,18,22,24]
        df = df[df['interact_topic'].isin(selected_int_topics)]
        dfmeta = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')
        df = pd.merge(df,dfmeta[['cell','subtype','cell_topic']],left_on=celltype,right_on='cell',how='left')

        df = df[df['subtype']==subtype]

        dfs = df.groupby(['cell_topic','interact_topic'])[celltype].count().reset_index()

        df_select_topic  = pd.DataFrame(dfmeta[dfmeta['cluster_celltype'].str.contains(celltype)]['cell_topic'].value_counts())
        selected_cell_topic=list(df_select_topic[df_select_topic['cell_topic']>100].index)

        dfs = dfs[dfs['cell_topic'].isin(selected_cell_topic)]

        celltopic_sum = dict(dfs.groupby('cell_topic')[celltype].sum())
        dfs['ncount'] = [x/celltopic_sum[y] for x,y in zip(dfs[celltype],dfs['cell_topic'])]

        dfs = dfs.rename(columns={celltype:'celltype_val'})
        dfs['celltype'] = subtype

        df_combined = pd.concat([df_combined, dfs], axis=0, ignore_index=True)

    df_combined.to_csv(spr.interaction_topic.id+'5_cancer_stype_nbr_ct_it_summary.csv.gz',index=False,compression='gzip')

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


    # first get neighbour cells of each cancer cell and their interaction topic 
    '''
    total data size will be 4185039
    total cancer cells 26321
            26321 x 159 neighbours = 4185039
    '''

    df=pd.read_csv(spr.interaction_topic.id+'3_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')

    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]
    dfmeta = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')
    df = pd.merge(df,dfmeta[['cell','cluster_celltype']],left_on='nbr',right_on='cell',how='left')

    # df = df[~df['cluster_celltype'].str.contains('Cancer')]

    experiment_config = read_config(experiment_home+'config.yaml')
    args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())

    spr.data.raw_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_data)


    corr_th = 0.3
    zcutoff = 4.0
    selected_int_topics = [2,4,7,10,18,22,24]
    

    ### correlation network

    dfall = _network.lr_correlation_network(spr,df,selected_int_topics,zcutoff,corr_th,'lr')
    dfall.to_csv(spr.interaction_topic.id+'8_lrnetwork.csv.gz',index=False,compression='gzip')


    ## chord network

    df_chord = _network.lr_chord_network(spr,df,selected_int_topics,zcutoff,corr_th,'lr')
    # df_chord = df_chord.sort_values(['ligand','receptor'])

    # df_chord['ligand'] = [x+'/'+str(y) for x,y in zip(df_chord['ligand'],df_chord['topic'])]
    # df_chord['receptor'] = [x+'/'+str(y) for x,y in zip(df_chord['receptor'],df_chord['topic'])]

    # df_chord['lrpair'] = [x+'_'+str(y) for x,y in zip(df_chord['ligand'],df_chord['receptor'])]

    # df_chord = df_chord.drop(columns=['ligand','receptor'])
    # df_chord = df_chord[['lrpair','topic']]

    # df_chordm = pd.melt(df_chord,id_vars=['topic'],value_vars=['ligand','receptor'])
    # df_chord['type']=['lr/'+str(x) for x in df_chord['topic']]
    # df_chordm['score'] = 1
    # df_chordm = df_chordm[['value','topic','variable']]


    # get and mark top genes
    from analysis import _topics
    top_n = 50
    df_top_genes = _topics.topic_top_lr_genes(spr,top_n)
    top_receptors = df_top_genes[df_top_genes['GeneType']=='ligands']['Gene'].unique()
    top_ligands = df_top_genes[df_top_genes['GeneType']=='receptors']['Gene'].unique()

    df_chord['ligand'] = [ x+'*' if x in top_ligands else x for x in df_chord['ligand']]
    df_chord['receptor'] = [ x+'*' if x in top_receptors else x for x in df_chord['receptor']]
    df_chord = df_chord.drop(columns=['score'])
    df_chord.to_csv(spr.interaction_topic.id+'8_chorddata.csv.gz',index=False,compression='gzip')

def correlation_lr_network_db():
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


    # first get neighbour cells of each cancer cell and their interaction topic 
    '''
    total data size will be 4185039
    total cancer cells 26321
            26321 x 159 neighbours = 4185039
    '''

    df=pd.read_csv(spr.interaction_topic.id+'3_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')

    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]
    dfmeta = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')
    df = pd.merge(df,dfmeta[['cell','cluster_celltype']],left_on='nbr',right_on='cell',how='left')

    # df = df[~df['cluster_celltype'].str.contains('Cancer')]

    experiment_config = read_config(experiment_home+'config.yaml')
    args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())

    spr.data.raw_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_data)


    corr_th = 0.3
    zcutoff = 4.0
    selected_int_topics = [2,4,7,10,18,22,24]
    

    ### correlation network

    df_db = pd.read_csv(experiment_home+args.database+args.lr_db,sep='\t')
    
    dfall = _network.lr_correlation_network_db(spr,df,selected_int_topics,zcutoff,corr_th,df_db)
    dfall.to_csv(spr.interaction_topic.id+'8_lrnetwork_lrdb.csv.gz',index=False,compression='gzip')

#################

def subtype_circle_plot():
    celltype='Cancer'
    df=pd.read_csv(spr.interaction_topic.id+'2_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')
    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]
    dfmeta = pd.read_csv(spr.interaction_topic.id+'2_a_meta_ct_argmax_maxprop.csv.gz')
    df = pd.merge(df,dfmeta,left_on=celltype,right_on='cell',how='left')

    # df = df[['Cancer','subtype','cluster_celltype','cell_topic','interact_topic_x']]


    st='TNBC'
    ct=7
    it=18

    st='HER2+'
    ct=24
    it=10

    nbr_cells = df[(df.subtype==st) & (df.cell_topic==ct) & (df.interact_topic_x ==it)].nbr.values
    dfmeta[dfmeta['cell'].isin(nbr_cells)].cluster_celltype.value_counts()

def cancer_nbr_normal_nonnormal():

    df = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')

    '''here we are using celltype instead of
    cluster_celltype because cancer cells are 
    present in non-normal cells only'''
    df = df[~df['celltype'].str.contains('Cancer')]

    # cancer vs normal in cell topic
    celltype='cell'
    dfs = df.groupby(['cell_topic','celltype'])[celltype].count().reset_index()

    dfs['cell_status'] = ['Normal' if '_n' in x else 'Cancer'for x in dfs['celltype']]

    dfs = dfs.groupby(['cell_topic','cell_status'])[celltype].sum().reset_index()

    celltopic_sum = dict(dfs.groupby('cell_topic')[celltype].sum())
    dfs['ncount'] = [x/celltopic_sum[y] for x,y in zip(dfs[celltype],dfs['cell_topic'])]

    dfs = dfs.rename(columns={celltype:'celltype_val'})

    # dfs = dfs[dfs['celltype_val']]
    dfs.to_csv(spr.interaction_topic.id+'4_cell_nbr_ct_norm_cancer_summary.csv.gz',index=False,compression='gzip')

    # cancer vs normal in interaction topic

    df=pd.read_csv(spr.interaction_topic.id+'3_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')
    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]

    celltype='Cancer'
    dfmeta = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')
    df = pd.merge(df,dfmeta,left_on='nbr',right_on='cell',how='left')

    '''here we are using celltype instead of
    cluster_celltype because cancer cells are 
    present in non-normal cells only'''
    df = df[~df['celltype'].str.contains('Cancer')]

    dfs = df.groupby(['interact_topic_x','celltype'])[celltype].count().reset_index()

    dfs['cell_status'] = ['Normal' if '_n' in x else 'Cancer'for x in dfs['celltype']]

    dfs = dfs.groupby(['interact_topic_x','cell_status'])[celltype].sum().reset_index()

    celltopic_sum = dict(dfs.groupby('interact_topic_x')[celltype].sum())
    dfs['ncount'] = [x/celltopic_sum[y] for x,y in zip(dfs[celltype],dfs['interact_topic_x'])]

    dfs = dfs.rename(columns={celltype:'celltype_val'})

    dfs.to_csv(spr.interaction_topic.id+'4_cell_nbr_it_norm_cancer_summary.csv.gz',index=False,compression='gzip')

def deg_analysis():
    ##################################################################
    # DEG analysis
    ##################################################################

    selected_int_topics = [2,4,10,18,22,24]

    df = pd.read_csv(spr.interaction_topic.id+'_it_model_all_topics.csv.gz',compression='gzip')
    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]

    df_its = pd.read_csv(spr.interaction_topic.id+'1_h_argmax.csv.gz',compression='gzip')
    df_nbr = spr.cell_topic.neighbour

    res = _topics.get_cell_neighbours_states(df_nbr,df_its,df['cell'].values)
    df_res = pd.DataFrame(res,columns=['cell_id','nbr','interact_topic'])
    df_res = df_res.explode(['nbr','interact_topic'])
    df_res = pd.merge(df_res,df[['cell','cell_topic','cluster','cluster_celltype']],left_on='nbr',right_on='cell',how='left')
    df_res = df_res.drop(columns=['cell'])

    df_res = df_res[df_res['interact_topic'].isin(selected_int_topics)]

    df_res.to_pickle(spr.interaction_topic.id+'_it_model_all_topics_all_cells.pkl')

    df_res = pd.read_pickle(spr.interaction_topic.id+'_it_model_all_topics_all_cells.pkl')

    query_cells = df_res['cell_id'].values
    query_cells = pd.Series(query_cells).unique()
    df = spr.interaction_topic_prop_with_cellids_nbrsummed(query_cells)
    df.to_csv(spr.interaction_topic.id+'_it_model_deg_analysis.csv.gz',index=False,compression='gzip')

    #### DEG data check
    df = pd.read_csv(spr.interaction_topic.id+'_it_model_deg_analysis.csv.gz')
    df_meta = pd.read_csv(spr.interaction_topic.id+'_alltopics_meta.csv.gz',compression='gzip')
    dftest = df_meta[df_meta['cell'].isin(df['cell_id'].values)]
    dftest.shape
    dftest['interact_topic'].value_counts()
    dftest.cluster_celltype.value_counts()

def heterogeneity_tnbc():
    
    celltype='Cancer'

    # print(subtype)
    df=pd.read_csv(spr.interaction_topic.id+'3_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')
    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]
    dfmeta = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')
    df = pd.merge(df,dfmeta[['cell','subtype','cell_topic']],left_on=celltype,right_on='cell',how='left')

    from util._io import read_config
    from collections import namedtuple
    experiment_config = read_config(experiment_home+'config.yaml')
    args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())
    dfraw = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_data)
    from anndata import AnnData
    import scanpy as sc
    import numpy as np
    adata = AnnData(dfraw.iloc[:,1:].values)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    dfn = adata.to_df()
    dfn.columns = dfraw.columns[1:]
    dfn['index'] = dfraw['index'].values
    dfn = dfn.set_index('index')

    ## select only cancer cells 
    dfn = dfn[dfn.index.isin(df['Cancer'].unique())]

    df = df.drop_duplicates(['Cancer','subtype','cell_topic'])

    dfjoin = pd.merge(dfn,df[['Cancer','subtype','cell_topic']],left_on=dfn.index,right_on='Cancer',how='left')

    dfjoin = dfjoin[ (     ( (dfjoin['cell_topic']==24) & (dfjoin['subtype'] =='HER2+') ) |    
                 ( (dfjoin['cell_topic']==48) & (dfjoin['subtype'] =='ER+')) |
                 ( (dfjoin['cell_topic']==19) & (dfjoin['subtype'] =='ER+')) |
                 ( (dfjoin['cell_topic']==9) & (dfjoin['subtype'] =='TNBC')) |
                 ( (dfjoin['cell_topic']==21) & (dfjoin['subtype'] =='TNBC')) |
                 ( (dfjoin['cell_topic']==24) & (dfjoin['subtype'] =='TNBC')) |
                 ( (dfjoin['cell_topic']==33) & (dfjoin['subtype'] =='TNBC')) |
                 ( (dfjoin['cell_topic']==43) & (dfjoin['subtype'] =='TNBC') )
            ) ]

    selected_int_topics = [24,22,18,10,7,4,2]
    tr,tl= _topics.get_topic_top_genes(spr,selected_int_topics,25)
    # dfn_l = dfn[tl]

    ext = ['Cancer','subtype','cell_topic']
    dfjoin_r = dfjoin[ext+tr]
    dfjoin_l = dfjoin[ext+tl]

    tr = list(pd.DataFrame(dfjoin_r.iloc[:,3:].mean(),columns=['val']).query("val>0.05").index)
    tl = list(pd.DataFrame(dfjoin_l.iloc[:,3:].mean(),columns=['val']).query("val>0.05").index)

    dfjoin_r = dfjoin[ext+tr]
    dfjoin_l = dfjoin[ext+tl]

    # dfjoin_l = dfjoin_l.groupby(['subtype','cell_topic']).sample(frac=0.1)
    # dfjoin_r = dfjoin_r.groupby(['subtype','cell_topic']).sample(frac=0.1)
    
    
    ### select topics for plot
    dfjoin_r.to_csv(spr.interaction_topic.id+'5_subtype_heatmap_r.csv.gz',index=False,compression='gzip')
    dfjoin_l.to_csv(spr.interaction_topic.id+'5_subtype_heatmap_l.csv.gz',index=False,compression='gzip')

def gse():
    from analysis import _gsea
    dfge = _gsea.gse_interactiontopic_lr_ranked(spr)
    # dfge.to_csv(spr.interaction_topic.id+'12_gseapy_out.csv',index=False)
    df.to_csv(spr.interaction_topic.id+'12_gsea.csv',index=False)
    
def survival():

    from analysis import _survival

    data='/data/TNBC/bulk/BRCA-US/'
    
    dfm = pd.read_csv(data+'donor.BRCA-US.tsv.gz',sep='\t')
    # dfm['overall_time'] = [ x if pd.isna(y) else y for x,y in zip(dfm['donor_interval_of_last_followup'],dfm['donor_survival_time'])]
    dfm = dfm[dfm['donor_vital_status']=='deceased']
    dfm['overall_time'] = dfm['donor_survival_time']
     
    df = pd.read_csv(data+'exp_seq.BRCA-US.tsv.gz',sep='\t')
    df = df[['icgc_donor_id','gene_id','normalized_read_count']]
    df = df[df['icgc_donor_id'].isin(dfm['icgc_donor_id'])]
    df = df.pivot_table(index=['icgc_donor_id'],columns='gene_id',values='normalized_read_count')

    df_score = _survival.generate_data_it(spr,df,dfm)
    df_score.to_csv(spr.interaction_topic.id+'11_survival_analysis_it_2g.csv.gz',index=False,compression='gzip')

    df_exp = np.exp(df_score.iloc[:,3:])
    df_exp = df_exp.div(df_exp.sum(axis=1), axis=0) 
    df_score[[2,4,7,10,18,22,24]] = df_exp
    df_score.to_csv(spr.interaction_topic.id+'11_survival_analysis_it_mod_1.csv.gz',index=False,compression='gzip')


    df_score = pd.melt(df_score,id_vars=['icgc_donor_id','donor_vital_status','overall_time'])                                            
    df_score.to_csv(spr.interaction_topic.id+'11_survival_analysis_it_mod.csv.gz',index=False,compression='gzip')
    
    
    df_score = _survival.generate_data_ct(spr,df,dfm)
    df_score.to_csv(spr.interaction_topic.id+'11_survival_analysis_ct.csv.gz',index=False,compression='gzip')

def supplemental_paper_revision():

    ## count total cell pair numbers 

    # df_kmeans = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')
    # df_its = pd.read_csv(spr.interaction_topic.id+'1_h_argmax.csv.gz',compression='gzip')
    # dfgrp = df_kmeans.groupby(['cluster_celltype','interact_topic']).count()['cell'].reset_index()
    # dfgrp = dfgrp.sort_values('cell',ascending=False)
    # dfgrp = dfgrp.drop_duplicates(['cluster_celltype','interact_topic'])
    # dfgrp = dfgrp.drop_duplicates(['interact_topic'])

    ## find 3 top pairs for each interaction topic

    ## get closest neighbours
    import csv
    df_its = pd.read_csv(spr.interaction_topic.id+'1_h_argmax.csv.gz',compression='gzip')
    df_its = df_its.iloc[:,:2]
    df_nbr = spr.cell_topic.neighbour.copy()

    top_pair = {}
    for idx,row in df_its.iterrows():
        it = row[1]
        cell = row[0]
        nbr = df_its.cell.values[df_nbr.iloc[idx,1]]
        if it not in top_pair.keys():
            top_pair[it] = [cell+'/'+nbr]
        else:
            if len(top_pair[it])<3:
                top_pair[it].append(cell+'/'+nbr)

    with open('dict.csv', 'w') as csv_file:  
        writer = csv.writer(csv_file)
        for key, value in top_pair.items():
            writer.writerow([key, value])


def paper_revision_correlation_lr_network_string_db():
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


    # first get neighbour cells of each cancer cell and their interaction topic 
    '''
    total data size will be 4185039
    total cancer cells 26321
            26321 x 159 neighbours = 4185039
    '''

    df=pd.read_csv(spr.interaction_topic.id+'3_meta_Cancer_cells_nbrs.csv.gz',compression='gzip')

    selected_int_topics = [2,4,7,10,18,22,24]
    df = df[df['interact_topic'].isin(selected_int_topics)]
    dfmeta = pd.read_csv(spr.interaction_topic.id+'3_a_meta_ct_argmax_maxprop.csv.gz')
    df = pd.merge(df,dfmeta[['cell','cluster_celltype']],left_on='nbr',right_on='cell',how='left')

    # df = df[~df['cluster_celltype'].str.contains('Cancer')]

    experiment_config = read_config(experiment_home+'config.yaml')
    args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())

    spr.data.raw_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_data)


    corr_th = 0.3
    zcutoff = 4.0
    selected_int_topics = [2,4,7,10,18,22,24]
    

    ### correlation network
    df_lrdb = pd.read_csv(experiment_home+args.database+args.lr_db,sep='\t')
    all_ligand = list(df_lrdb['ligand_gene_symbol'].unique())
    all_receptor = list(df_lrdb['receptor_gene_symbol'].unique())

    df_db = pd.read_csv(experiment_home+args.data+'string_db/string_db.csv.gz')
    df_db = df_db[df_db['score']>800]
    df_db = df_db[df_db['gene1'].isin(all_ligand)]
    df_db = df_db[df_db['gene2'].isin(all_receptor)]
    df_db.columns = ['ligands','receptors','score']


    dfall = _network.lr_correlation_network_stringdb(spr,df,selected_int_topics,zcutoff,corr_th,df_db)
    dfall.to_csv(spr.interaction_topic.id+'8_lrnetwork_stringdb.csv.gz',index=False,compression='gzip')
