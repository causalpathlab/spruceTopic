import experiment_cell_topic 
import pandas as pd

experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'
spr = experiment_cell_topic.get_experiment(experiment_home)
print(spr.cell_topic.model_id)

def generate_umap_coordinates():

    # 1.1 generate and save umap coordinates based on softmaxed latent dimensions 

    import umap

    dfh = spr.cell_topic.h.copy()
    df_umap= pd.DataFrame()
    df_umap['cell'] = dfh['cell']

    umap_2d = umap.UMAP(n_components=2, init='random', random_state=0,min_dist=0.0,metric='cosine')
    proj_2d = umap_2d.fit(dfh.iloc[:,1:])

    df_umap[['umap1','umap2']] = umap_2d.embedding_[:,[0,1]]
    df_umap['cell_topic'] = [x.replace('h','') for x in dfh.iloc[:,1:].idxmax(axis=1)]
    df_umap.to_csv(spr.cell_topic.id+'1_a_umap_cordinates.csv.gz',index=False,compression='gzip')

def umap_plots_argmax():

    ##################################################################
    # 1.2 
    # get previously saved umap coordinates and
    # generate argmax label plot

    import matplotlib.pylab as plt
    plt.rcParams['figure.figsize'] = [15, 10]
    plt.rcParams['figure.autolayout'] = True
    import colorcet as cc
    import seaborn as sns

    df_umap = pd.read_csv(spr.cell_topic.id+'1_a_umap_cordinates.csv.gz',compression='gzip')

    # if we need to select topics with >100 cells
    # selected_topic = list(df_umap.topic.value_counts().index[0:32])
    # df_umap = df_umap[df_umap['cell_topic'].isin(selected_topic)]

    # plot umap with argmax topic

    cp = sns.color_palette(cc.glasbey_dark, n_colors=len(df_umap['cell_topic'].unique()))
    p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='cell_topic',s=2,palette=cp,legend=False)
    # plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
    p.axes.set_title("UMAP plot of cell-mixture with argmax cell topic assignment",fontsize=30)
    p.set_xlabel("UMAP1",fontsize=20)
    p.set_ylabel("UMAP2",fontsize=20)
    plt.savefig(spr.cell_topic.id+'1_a_umap_ct_argmax.png',dpi=600);plt.close()

def kmeans_label():

    from sklearn.cluster import KMeans

    df_umap = pd.read_csv(spr.cell_topic.id +'1_a_umap_cordinates.csv.gz')

    dfh = spr.cell_topic.h.copy()

    kmeans = KMeans(n_clusters=50, random_state=0).fit(dfh.iloc[:,1:].to_numpy())

    df_kmeans = dfh[['cell']]
    df_kmeans['cluster'] = kmeans.labels_
    df_umap = pd.merge(df_umap,df_kmeans,on='cell',how='left')

    ##### save 
    df_umap.to_csv(spr.cell_topic.id+'1_c_kmeans.csv.gz',index=False,compression='gzip')

def umap_add_cell_annotation_label():

    ##################################################################
    # 1.3 
    # get previously saved umap coordinates and
    # generate umap plot with celltype label from metadata and annotated files

    df_umap = pd.read_csv(spr.cell_topic.id+'1_c_kmeans.csv.gz',compression='gzip')
    df_umap['label'] = [x.split('_')[len(x.split('_'))-1] for x in df_umap['cell']]

    # metadata
    f='/home/BCCRC.CA/ssubedi/projects/data/GSE176078mix/GSE176078_metadata.csv.gz'
    dfl = pd.read_csv(f,compression='gzip')
    dfl = dfl.rename(columns={'Unnamed: 0':'cell'})

    dflabel = pd.DataFrame()
    dflabel['l1'] =  [x for x in df_umap[df_umap['label']=='GSE176078']['cell']]
    dflabel['l2'] =  [x.replace('_GSE176078','') for x in df_umap[df_umap['label']=='GSE176078']['cell']]
    dflabel = pd.merge(dflabel,dfl,right_on='cell',left_on='l2',how='left')

    print(dfl.columns)
    labels=['l1','subtype', 'celltype_subset', 'celltype_minor', 'celltype_major']

    df_umap = pd.merge(df_umap,dflabel[labels],right_on='l1',left_on='cell',how='left')
    for label in labels:df_umap[label] = df_umap[label].mask(df_umap[label].isna(), df_umap['label'])


    label='celltype_major'
    # annotated file
    normal_f='/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/output/scanpy/gse_normal/gse_normal_annotated_chetah.csv.gz'
    df_nlabel =  pd.read_csv(normal_f,compression='gzip')
    df_nlabel = df_nlabel[~df_nlabel['pruned_labels'].isna()]
    df_nlabel['pruned_labels'] = [ x.replace(' ','_')+'_n' for x in df_nlabel['pruned_labels']]

    df_umap = pd.merge(df_umap,df_nlabel[['cell','pruned_labels']],on='cell',how='left')
    df_umap.loc[df_umap['pruned_labels'].notna(),label] = df_umap[df_umap['pruned_labels'].notna()]['pruned_labels'].values
    df_umap =df_umap.drop(columns=['pruned_labels'])
    df_umap.loc[df_umap[label]=='GSE164898',label] = 'Epithelial_n'
    df_umap.loc[df_umap[label]=='GSE156728-CD8',label] = 'T_CD8_pan'
    df_umap.loc[df_umap[label]=='GSE156728-CD4',label] = 'T_CD4_pan'
    df_umap[label] = df_umap[label].mask(df_umap[label].isna(), df_umap['label'])

    # fix label names for visualization

    df_umap[label] = [ x.replace('-','_').replace(' ','_')  for x in df_umap[label]]

    df_umap = df_umap.sort_values(label)


    df_clust_celltype = df_umap.groupby(['cluster',label])['cell'].count().reset_index().sort_values(['cluster','cell'])
    df_clust_celltype = df_clust_celltype.drop_duplicates(subset=['cluster'], keep='last').drop(['cell'],axis=1).rename(columns={label:'cluster_celltype'})

    df_umap = pd.merge(df_umap,df_clust_celltype,on='cluster',how='left')

    clust_celltype_group={
    'B_cell_n':'B',
    'B_cells':'B',

    'CAF_n':'CAF',
    'CAFs':'CAF',
    'Myofibroblast_n':'CAF',

    'Cancer_Epithelial':'Cancer',

    'Endothelial':'Endothelial',
    'Endothelial_n':'Endothelial',

    'Epithelial_n':'Epithelial',
    'Normal_Epithelial':'Epithelial',

    'Macrophage_n':'Myeloid',
    'Myeloid':'Myeloid',
    'Monocyte':'Myeloid',

    'PVL':'PVL',

    'Plasma_n':'Plasma',
    'Plasmablasts':'Plasma',

    'T_CD8_pan':'T',
    'T_cells':'T',
    }


    df_umap['celltype'] = [ clust_celltype_group[x] for x in df_umap['cluster_celltype']]

    df_umap['celltype_ct'] = [x+'/'+str(y) for x,y in zip(df_umap['celltype'],df_umap['cell_topic'])]

    df_umap = df_umap.sort_values('celltype')

    df_umap.to_csv(spr.cell_topic.id +'1_d_celltopic_label.csv.gz',index=False,compression='gzip')

#### plot figures ####

def umap_plots_kmeans_label():

    import matplotlib.pylab as plt
    plt.rcParams['figure.figsize'] = [15, 10]
    plt.rcParams['figure.autolayout'] = True
    import colorcet as cc
    import seaborn as sns

    df_umap = pd.read_csv(spr.cell_topic.id +'1_d_celltopic_label.csv.gz')

    df_umap = df_umap.sort_values('celltype')

    colors =  ["orange", "#1CE6FF", "fuchsia", "#7A4900", "green","#006FA6" ,"limegreen", "red", "blue"]

    cp = sns.color_palette(colors,9)
    p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='celltype',s=2,palette=cp,legend=False)
    # plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
    # p.axes.set_title("UMAP plot of cell-mixture kmeans label",fontsize=30)
    p.set_xlabel("UMAP1",fontsize=20)
    p.set_ylabel("UMAP2",fontsize=20)
    plt.axis('off')
    plt.savefig(spr.cell_topic.id+'2_a_kmeans_celltype_grouped.png',dpi=600);plt.close()

def latent_hist_plots():
    ##################################################################
    # 1 a 
    ##################################################################
    import matplotlib.pylab as plt
    plt.rcParams['figure.figsize'] = [20, 5]
    plt.rcParams['figure.autolayout'] = True
    import seaborn as sns

    dfh = spr.cell_topic.h.copy()
    dfh = dfh.set_index('cell')
    dfhm = dfh.melt()
    dfhm['variable'] = [x.replace('h','') for x in dfhm['variable']]
    dfhm.columns = ['Cell topic','Topic proportion distribution']
    p = sns.stripplot(x='Cell topic', y='Topic proportion distribution', data=dfhm,jitter=True, linewidth=0.1,size=0.5)
    _ = plt.xticks(rotation=45, ha='right')
    p.set_xlabel('Cell topic',fontsize=20)
    p.set_ylabel('Topic proportion distribution',fontsize=20)
    plt.savefig(spr.cell_topic.id+'2_a_topic_distribution.png',dpi=600);plt.close()

    ##################################################################
    # 1 b
    ##################################################################
    import matplotlib.pylab as plt
    plt.rcParams['figure.figsize'] = [20, 5]
    plt.rcParams['figure.autolayout'] = True
    import seaborn as sns

    dfh = spr.cell_topic.h.copy()
    df_hmax = pd.DataFrame(pd.Series([x.replace('h','') for x in dfh.iloc[:,1:].idxmax(axis=1)]).value_counts()).reset_index().rename(columns={'index':'cell_topic',0:'argmax_count'})
    df_hmax = df_hmax.sort_values('cell_topic')
    p = sns.barplot(x='cell_topic',y='argmax_count',data=df_hmax,palette=sns.color_palette("hls", 50))
    p.set_xlabel("Cell Topic",fontsize=20)
    p.set_ylabel("Total cells",fontsize=20)
    plt.savefig(spr.cell_topic.id+'2_a_ct_argmax.png',dpi=600);plt.close()

    '''    In [73]: df_hmax.argmax_count.mean()
        Out[73]: 3118.26

        In [74]: df_hmax.argmax_count.std()
        Out[74]: 5987.42921157933
    
    We found each cell topic corresponds to a group of average 3118 cells (with standard deviation $\pm$ 5987).
    
    '''

def latent_hist_plots_withcelltype():

    #### plot major cell types

    import matplotlib.pylab as plt
    plt.rcParams['figure.figsize'] = [20, 5]
    plt.rcParams['figure.autolayout'] = True
    import colorcet as cc
    import seaborn as sns
    sns.set_style("white",{'axes.linewidth': 1})
    import numpy as np

    df = pd.read_csv(spr.cell_topic.id +'1_d_celltopic_label.csv.gz')


    dfs = df.groupby(['cell_topic','celltype'])['cell'].count().reset_index()
    celltopic_sum = dict(dfs.groupby('cell_topic')['cell'].sum())
    dfs['ncount'] = [x/celltopic_sum[y] for x,y in zip(dfs['cell'],dfs['cell_topic'])]

    dfs = dfs.rename(columns={'cell':'celltype_val'})
    dfs = dfs.sort_values('celltype')

    colors =  ["orange", "#1CE6FF", "fuchsia", "#7A4900", "green","#006FA6" ,"limegreen", "red", "blue"]

    ax=sns.histplot(data=dfs,x="cell_topic",
    hue="celltype",weights="ncount",multiple="stack",discrete=True,palette=sns.color_palette(colors, 9))
    ax.set(xticks=np.arange(0,50,1))
    ax.set_xlabel("Cell topic",fontsize=20)
    ax.set_ylabel("Cell type distribution",fontsize=20)
    sns.move_legend(ax, "center left", bbox_to_anchor=(1, 0.5))
    plt.savefig(spr.cell_topic.id+'2_b_ct_celltype_dist_norm.pdf',format="pdf",dpi=600);plt.close()

    ax=sns.histplot(data=dfs,x="cell_topic",
    hue="celltype",weights="celltype_val",multiple="stack",discrete=True,palette=sns.color_palette(colors, 9))
    ax.set(xticks=np.arange(0,50,1))
    ax.set_xlabel("Cell topic",fontsize=20)
    ax.set_ylabel("Cell type distribution",fontsize=20)
    sns.move_legend(ax, "center left", bbox_to_anchor=(1, 0.5))
    plt.savefig(spr.cell_topic.id+'2_c_ct_celltype_dist.pdf',format="pdf",dpi=600);plt.close()

    dfs.to_csv(spr.cell_topic.id+'2_d_ct_celltype_dist_norm.csv',index=False)

def umap_plots_kmeans_label():

    ##################################################################
    # 1.2 
    # get previously saved umap coordinates and
    # generate argmax label plot

    import matplotlib.pylab as plt
    plt.rcParams['figure.figsize'] = [15, 10]
    plt.rcParams['figure.autolayout'] = True
    import colorcet as cc
    import seaborn as sns

    df_umap = pd.read_csv(spr.cell_topic.id+'1_d_celltopic_label.csv.gz',compression='gzip')

    # if we need to select topics with >100 cells
    # selected_topic = list(df_umap.topic.value_counts().index[0:32])
    # df_umap = df_umap[df_umap['cell_topic'].isin(selected_topic)]

    # plot umap with argmax topic

    cp = sns.color_palette(cc.glasbey_dark, n_colors=len(df_umap['cell_topic'].unique()))
    p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='cell_topic',s=2,palette=cp,legend=False)
    # plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
    p.axes.set_title("UMAP plot of cell-mixture with argmax cell topic assignment",fontsize=30)
    p.set_xlabel("UMAP1",fontsize=20)
    p.set_ylabel("UMAP2",fontsize=20)
    plt.savefig(spr.cell_topic.id+'2_c_umap_ct_argmax.png',dpi=600);plt.close()


    df_umap = pd.read_csv(spr.cell_topic.id +'1_d_celltopic_label.csv.gz')
    # plot all labels
    import matplotlib.pylab as plt
    plt.rcParams['figure.figsize'] = [15, 10]
    plt.rcParams['figure.autolayout'] = True
    import colorcet as cc
    import seaborn as sns

    cp = sns.color_palette(cc.glasbey, n_colors=len(df_umap['celltype'].unique()))
    p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='celltype',s=2,palette=cp)
    plt.legend(title='Cell type',title_fontsize=18, fontsize=8,loc='center left', bbox_to_anchor=(1, 0.5))
    p.axes.set_title("UMAP plot of cell-mixture with cell type label",fontsize=30)
    # plt.figtext(.1,.8,"singleR/CHETAH for normal dataset",fontsize=20)
    p.set_xlabel("UMAP1",fontsize=20)
    p.set_ylabel("UMAP2",fontsize=20)
    plt.savefig(spr.cell_topic.id+'2_d_umap_celltype.png',dpi=600);plt.close()

    #### plot major cell types

    import matplotlib.pylab as plt
    plt.rcParams['figure.figsize'] = [15, 10]
    plt.rcParams['figure.autolayout'] = True
    import colorcet as cc
    import seaborn as sns

    df_umap = pd.read_csv(spr.cell_topic.id +'1_d_celltopic_label.csv.gz')

    df_umap = df_umap.sort_values('celltype')
    # cp = sns.color_palette(cc.glasbey, n_colors=len(df_umap['celltype'].unique()))

    colors =  ["orange", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941","#006FA6","yellowgreen" ,"#A30059","#FFDBE5", "#7A4900", "#0000A6"]
    sns.set_palette(sns.color_palette(colors))
    p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='celltype',s=2,palette=sns.color_palette(colors, 11))
    plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
    p.axes.set_title("UMAP plot of cell-mixture with major cell type label",fontsize=30)
    # plt.figtext(.1,.8,"singleR/CHETAH for normal dataset",fontsize=20)
    p.set_xlabel("UMAP1",fontsize=20)
    p.set_ylabel("UMAP2",fontsize=20)
    plt.savefig(spr.cell_topic.id+'2_e_umap_celltype_grouped.png',format='png',dpi=600);plt.close()

def umap_plot_cell_subtype_label():

    df_umap = pd.read_csv(spr.cell_topic.id +'1_d_celltopic_label.csv.gz')
    df_umap['selected'] = [ x  if x in ['TNBC','ER+','HER2+'] else 'Others' for x in df_umap['subtype']]
    df_umap = df_umap.sort_values('selected')

    colors = ["Red","Orange","Gray", "Green"]
    sns.set_palette(sns.color_palette(colors))
    p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='selected',s=0.2)
    plt.legend(title='Subtype',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
    p.axes.set_title("UMAP plot of cell-mixture with subtype label",fontsize=30)
    # plt.figtext(.1,.8,"singleR/CHETAH for normal dataset",fontsize=20)
    p.set_xlabel("UMAP1",fontsize=20)
    p.set_ylabel("UMAP2",fontsize=20)
    plt.savefig(spr.cell_topic.id+'2_f_umap_celltype_subtype.png',format='png',dpi=600);plt.close()

def umap_plot_cell_normal_dataset_label():

    df_umap = pd.read_csv(spr.cell_topic.id +'1_d_celltopic_label.csv.gz')
    df_umap['selected'] = [ x  if '_n' in x else 'Others' for x in df_umap['cluster_celltype']]
    df_umap = df_umap.sort_values('selected')
    colors =  ["orange", "#1CE6FF", "fuchsia", "#7A4900", "green","gray","#006FA6" ,"limegreen", "red", "blue"]
    colors =  ["orange", "#1CE6FF",  "#7A4900","green","#006FA6" ,"gray", "red",]


    sns.set_palette(sns.color_palette(colors))
    p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='selected',s=0.2)
    plt.legend(title='Celltype',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
    p.axes.set_title("UMAP plot of cells from normal dataset",fontsize=30)
    # plt.figtext(.1,.8,"singleR/CHETAH for normal dataset",fontsize=20)
    p.set_xlabel("UMAP1",fontsize=20)
    p.set_ylabel("UMAP2",fontsize=20)
    plt.savefig(spr.cell_topic.id+'2_g_umap_celltype_gsenormal.png',format='png',dpi=600);plt.close()


    df_umap['selected'] = [ x  if '_n' not in x else 'Others' for x in df_umap['cluster_celltype']]
    df_umap['selected'] = [ x  if '_pan' not in x else 'Others' for x in df_umap['selected']]
    df_umap = df_umap.sort_values('selected')
    colors =  ["orange", "#1CE6FF", "fuchsia", "#7A4900", "#006FA6" ,"green","gray","limegreen", "red", "blue"]


    sns.set_palette(sns.color_palette(colors))
    p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='selected',s=0.2)
    plt.legend(title='Celltype',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
    p.axes.set_title("UMAP plot of cells from cancer dataset",fontsize=30)
    # plt.figtext(.1,.8,"singleR/CHETAH for normal dataset",fontsize=20)
    p.set_xlabel("UMAP1",fontsize=20)
    p.set_ylabel("UMAP2",fontsize=20)
    plt.savefig(spr.cell_topic.id+'2_g_umap_celltype_gsecancer.png',format='png',dpi=600);plt.close()

def umap_plot_cell_pan_dataset_label():

    df_umap = pd.read_csv(spr.cell_topic.id +'1_d_celltopic_label.csv.gz')
    selected = [ x  if '_pan' in x else 'Others' for x in df_umap['celltype_major']]
    colors = ["Gray","Green", "Red"]
    sns.set_palette(sns.color_palette(colors))
    p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue=selected,s=0.2)
    plt.legend(title='Celltype',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
    p.axes.set_title("UMAP plot of cells from pan dataset",fontsize=30)
    # plt.figtext(.1,.8,"singleR/CHETAH for normal dataset",fontsize=20)
    p.set_xlabel("UMAP1",fontsize=20)
    p.set_ylabel("UMAP2",fontsize=20)
    plt.savefig(spr.cell_topic.id+'2_h_umap_celltype_pan.png',format="png",dpi=600);plt.close()

def umap_plot_cell_kmeans_cluster_label():

    df_umap = pd.read_csv(spr.cell_topic.id +'1_d_celltopic_label.csv.gz')
    ######

    cp = sns.color_palette(cc.glasbey, n_colors=len(df_umap['cluster'].unique()))
    p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='cluster',s=2,palette=cp,legend=False)

    # plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
    p.axes.set_title("UMAP plot of cell-mixture kmeans label",fontsize=30)
    p.set_xlabel("UMAP1",fontsize=20)
    p.set_ylabel("UMAP2",fontsize=20)
    plt.savefig(spr.cell_topic.id+'2_j_kmeans_cluster.png',dpi=600);plt.close()

def struct_plot_from_kmeans():
    ##################################################################
    # 2 data for struct plot in R
    ##################################################################

    dfh = spr.cell_topic.h.copy()
    df_umap = pd.read_csv(spr.cell_topic.id+'1_d_celltopic_label.csv.gz',compression='gzip')

    df_umap_sum = df_umap.groupby(['cell_topic','celltype_ct'])['cell'].count().reset_index()
    df_umap_sum = df_umap_sum.sort_values(['cell_topic','cell'])
    
    # important keep highest cell count to each topic
    df_umap_sum = df_umap_sum.drop_duplicates(subset=['cell_topic'], keep='last')

    df_umap_sum = df_umap_sum[df_umap_sum['cell']>100]

    df_umap_sum = df_umap_sum.sort_values('celltype_ct')

    selected_celltype_ct = df_umap_sum['celltype_ct'].values

    df_umap = df_umap[df_umap['celltype_ct'].isin(selected_celltype_ct)]

    df_h_sample_kmeans = df_umap.groupby('celltype_ct').sample(n=50, random_state=1)
    dfh = pd.merge(dfh,df_h_sample_kmeans[['cell','celltype_ct']],on='cell',how='inner')
    dfh = dfh.rename(columns={'celltype_ct':'Topic'})
    dfh.to_csv(spr.cell_topic.id+'2_i_kmeans_celltype_sample.csv.gz',index=False,compression='gzip')
    del dfh

    dfh = spr.cell_topic.h.copy()
    df_umap = pd.read_csv(spr.cell_topic.id+'2_i_kmeans.csv.gz',compression='gzip')

    df_h_sample_kmeans = df_umap.groupby('cluster').sample(n=50, random_state=1)
    dfh = pd.merge(dfh,df_h_sample_kmeans[['cell','cluster']],on='cell',how='inner')
    dfh = dfh.rename(columns={'cluster':'Topic'})
    dfh.to_csv(spr.cell_topic.id+'2_i_kmeans_cluster_sample.csv.gz',index=False,compression='gzip')
    del dfh

def plot_top_genes():
    ##################################################################
    # 4 data for top genes heatmap in R 
    ##################################################################
    
    from analysis import _topics

    top_n = 25
    df_top_genes = _topics.topic_top_genes(spr,top_n)
    df_top_genes.to_csv(spr.cell_topic.id+'4_beta_weight_top_'+str(top_n)+'_genes.csv.gz',index=False,compression='gzip')


    #######
    df_umap = pd.read_csv(spr.cell_topic.id +'1_d_celltopic_label.csv.gz')

    df_umap_sum = df_umap.groupby(['cell_topic','celltype_ct'])['cell'].count().reset_index()
    df_umap_sum = df_umap_sum.sort_values(['cell_topic','cell'])
    
    # important keep highest cell count to each topic
    df_umap_sum = df_umap_sum.drop_duplicates(subset=['cell_topic'], keep='last')

    df_umap_sum = df_umap_sum[df_umap_sum['cell']>100]

    df_umap_sum = df_umap_sum.sort_values('celltype_ct')


    top_n = 25
    df_top_genes = _topics.topic_top_genes(spr,top_n)
    sel_topics = ['k'+str(x) for x in df_umap_sum['cell_topic']]
    df_top_genes = df_top_genes[df_top_genes['Topic'].isin(sel_topics)]

    df_top_genes = df_top_genes.pivot(index='Topic',columns='Gene',values='Proportion')
    df_top_genes[df_top_genes>10] = 10
    df_top_genes[df_top_genes<-10] = -10

    topic_label = [ x for x in zip(df_umap_sum['cell_topic'],df_umap_sum['celltype_ct'])]
    tld = {}
    for x in topic_label:tld[x[0]]= x[1]

    df_top_genes.index = [ tld[int(x.replace('k',''))] for x in df_top_genes.index]
    
    df_top_genes.to_csv(spr.cell_topic.id+'4_beta_weight_top_'+str(top_n)+'_genes_selected.csv.gz',index=True,compression='gzip')

    top3=[]
    for i in range(df_top_genes.shape[0]):top3.append(df_top_genes.iloc[i,:].sort_values(ascending=False).index[:3].values)

    df_top3 = pd.DataFrame(top3,columns=['a','b','c'])
    df_top3['genes'] =  df_top3['a'] +'/' +df_top3['b']+'/'+df_top3['c']
    df_top3.index = df_top_genes.index
    df_top3 = df_top3.drop(columns=['a','b','c'])
    df_top3.to_csv(spr.cell_topic.id+'4_beta_weight_top_'+str(top_n)+'_genes_top3.csv.gz',index=True,compression='gzip')

def plot_marker_genes():
    #### plot marker genes

    from util._io import read_config
    from collections import namedtuple


    experiment_home='/home/BCCRC.CA/ssubedi/'+experiment_home
    experiment_config = read_config(experiment_home+'config.yaml')
    args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())
    df = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_data)

    from anndata import AnnData
    import scanpy as sc
    import numpy as np


    adata = AnnData(df.iloc[:,1:].values)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    dfn = adata.to_df()
    dfn.columns = df.columns[1:]
    dfn['cell'] = df['index'].values

    df_umap = pd.read_csv(spr.cell_topic.id+'1_a_umap_cordinates.csv.gz',compression='gzip')

    dfn = pd.merge(dfn,df_umap,on='cell',how='left')
    # df_m = df.groupby('topic').sample(frac=0.1, random_state=1)
    # dfn = dfn[dfn['cell'].str.contains('GSE176078')]

    from analysis import _supp
    # marker_genes = [ 'CD79A', 'MS4A1','CD3E', 'CD3G', 'CD3D','CD4', 'CD8A', 'CD8B',
    # 'SELL', 'CCR7', 'IL7R']
    # marker_genes = ['S100A8', 'S100A9', 'CST3','FCGR3A', 'FCGR3B', 'CD14'] 
    # marker_genes = ['EPCAM', 'MKI67', 'CD3D', 'CD68','MS4A1','JCHAIN', 'PECAM1', 'PDGFRB'] 
    marker_genes = ['EPCAM', 'CD3D', 'CD68','MS4A1', 'PECAM1', 'PDGFRB'] 

    _supp.plot_marker(spr,dfn,marker_genes)
