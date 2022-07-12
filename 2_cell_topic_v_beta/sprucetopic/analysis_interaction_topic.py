##################################################################
    # setup
##################################################################
import experiment_interaction_topic 
import pandas as pd
from pathlib import Path

server = Path.home().as_posix()
experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'
experiment_home = server+experiment_home
spr = experiment_interaction_topic.get_experiment_model(experiment_home)
print(spr.cell_topic.model_id)
print(spr.interaction_topic.model_id)

##################################################################
# 1 analysis of latent h
##################################################################
from analysis import _topics
import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [15, 10]
plt.rcParams['figure.autolayout'] = True
import colorcet as cc
import seaborn as sns

'''
The latent h of interaction topic is based on cell-pairs. So, we have to take argmax for each pair. 
Each row is one cell and argmax h for its 159 neighbours.
'''
df_its = spr.interaction_topic_states()

###save/load
df_its.to_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz',index=False,compression='gzip')
df_its = pd.read_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz',compression='gzip')

df_hmax = pd.DataFrame(pd.Series(df_its.iloc[:,1:].values.flatten()).value_counts()).reset_index().rename(columns={'index':'topic',0:'argmax_count'})
df_hmax = df_hmax.sort_values('argmax_count',ascending=False)
p = sns.barplot(x='topic',y='argmax_count',data=df_hmax,color='blue')
p.set_xlabel("Topic",fontsize=20)
p.set_ylabel("Count(argmax)",fontsize=20)
plt.savefig(spr.interaction_topic.model_id+'_it_h_argmax.png');plt.close()



##################################################################
# 2 top genes
##################################################################


from analysis import _topics

#### get top genes based on beta weight matrix - l/r separately

top_n = 10
df_top_genes = _topics.topic_top_lr_genes(spr,top_n)
df_top_genes.to_csv(spr.interaction_topic.model_id+'_it_beta_weight_top_'+str(top_n)+'_genes.csv.gz',index=False,compression='gzip')

#### get top genes based on beta weight matrix - l/r together

# df_db = pd.read_csv( experiment_home + spr.args.database+ spr.args.lr_db,sep='\t', usecols=['lr_pair'])
# top_n=25
# df_lr_topic = _topics.topic_top_lr_pair_genes(spr,df_db,top_n)
# df_lr_topic.to_csv(spr.interaction_topic.model_id+'_it_beta_top_'+str(top_n)+'_lrpair.csv.gz',compression='gzip')

##################################################################
# 8 interaction topic struct plot and cancer nbr gene expression
##################################################################

##################################################################
# 8.1 
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
df_kmeans = df_kmeans[ df_kmeans['topic'].isin(df_kmeans.topic.value_counts().index[:13])]
df_kmeans = df_kmeans.groupby('topic').sample(n=100, random_state=1)
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

iofit = [3,11,13,16,18,20,23]

dfh_sample_kmeans = dfh_sample_kmeans[dfh_sample_kmeans['Topic'].isin(iofit)]

dfh_sample_kmeans = dfh_sample_kmeans.sort_values('Topic')
dfh_sample_kmeans.to_csv(spr.interaction_topic.model_id+'_it_h_sample_kmeans_selected.csv.gz',index=False,compression='gzip')

##################################################################
# 8.2 
# get mean expression of cancer nbr cells from above sample dataset
# draw heatmap
'''
-for each interaction topic get cancer cells select its neighbour cells
-get raw data for these cancer cells and nbr cells
- normalize rowwise such that each ligand expression is proportion of all ligand expression
- get top 25 ligand and receptor genes for interaction topic
- filter those top genes and take mean, add 1e-5, and apply sqrt for heatmap
'''

df=pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics_cancer_cells.csv.gz',compression='gzip')
df = df[['cancer_cell','nbr','state']]
iofit = [2,4,7,10,18,22,24]
df = df[df['state'].isin(iofit)]
_topics.plt_cn(spr,df)


########### cell type distribution of cancer cells
# it with interesting pattern - 2,4,7,10,18,22,24


iofit = [2,4,7,10,18,22,24]
df=pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics_cancer_cells.csv.gz',compression='gzip')
df = df[~df['cluster_celltype'].str.contains('Cancer')]
df = df[df['state'].isin(iofit)]

df_grp = df.groupby(['cluster_celltype','state'])['state'].size().rename('count').reset_index()
##normalize
celltype_sum = dict(df_grp.groupby('state')['count'].sum())

df_grp['ncount'] = [x/celltype_sum[y] for x,y in zip(df_grp['count'],df_grp['state'])]

df_grp.to_csv(spr.interaction_topic.model_id+'_it_cancercells_interactions.csv.gz',index=False,compression='gzip')

##################################################################
# 3 cancer centric view
##################################################################

##################################################################
# 3.1 
# generate topics summary 

'''
Take umap coordinates, topic assignment, cluster assignment, cell type, and cluster majority cell type from
cell topic analysis and merge with state assignment from interaction topic.
'''
spr.interaction_topic.neighbour_h =  pd.read_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz')
df_kmeans = pd.read_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_kmeans.csv.gz')

df_all_topics = _topics.get_topics(spr,df_kmeans)
df_all_topics.to_csv(spr.interaction_topic.model_id+'_it_model_all_topics.csv.gz',index=False,compression='gzip')

##################################################################
# 3.2 
# look at all the cancer cells, what cell and interaction topics they belong to
# check only for interaction with >100 cells 

'''
Generate a model summary file where each cluster is assigned-
1. cell topic based on argmax
2. interaction state based on highest state assigned among all neighbours
3. total number of cells in the cluster
'''
df_kmeans = pd.read_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_kmeans.csv.gz')
df_kmeans = df_kmeans[ df_kmeans.cluster_celltype.str.contains('Cancer')]

dfsummary = _topics.topics_summary(spr,df_kmeans)
dfsummary = dfsummary.groupby(['cluster','topic','state']).count()
dfsummary = dfsummary.reset_index()	
### filter topics summary for summary plot

# use cell type assignment from label and annotation
clust_to_cell_type = dict(zip( df_kmeans['cluster'],df_kmeans['cluster_celltype']))
dfsummary['celltype'] = [clust_to_cell_type[x] for x in dfsummary['cluster']]

# use cluster assignment as cell type
# dfsummary['celltype'] = ['ct_'+str(x) for x in dfsummary['cluster']]

dfsummary = dfsummary.sort_values('state')
min_cells_per_state = 100
dfsummary = dfsummary[dfsummary['cell']>min_cells_per_state]
dfsummary.to_csv(spr.interaction_topic.model_id+'_it_model_summary.csv.gz',index=False,compression='gzip')



##################################################################
# 3.2 
# interaction topics distribution of neighbouring cells of cancer

# 3.2.1 get neighbour cells of each cancer cell and their interaction topic 

from analysis import _topics
df = pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics.csv.gz',compression='gzip')
df_its = pd.read_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz',compression='gzip')
df_nbr = spr.cell_topic.neighbour

res = _topics.get_cell_neighbours_states(df_nbr,df_its,df)
df_res = pd.DataFrame(res,columns=['cancer_cell','nbr','state'])
df_res = df_res.explode(['nbr','state'])
df_res = pd.merge(df_res,df[['cell','cluster','cluster_celltype']],left_on='nbr',right_on='cell',how='left')
df_res.to_csv(spr.interaction_topic.model_id+'_it_model_all_topics_cancer_cells.csv.gz',index=False,compression='gzip')


# 3.2.1 look at non cancer neighbours of cancer cells and check state distribution 

df = df_res


##################################################################
# 4 network graph
# Check interaction state of cancer cells and neighbouring cells 

'''
take all cancer cells and remove state 7 and 10 
plot network to show cancer cell states and neighbour cells belonging
to that state
'''
from analysis import _network
df = pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics_cancer_cells.csv.gz',compression='gzip')
# df = df[df['cluster_celltype'] == 'Cancer Epithelial']
# df = df[~df['state'].isin([7,10])]
_network.cell_interaction_network(spr,df)


##################################################################
# 5 topic wise lr network graph
##################################################################

from analysis import _network

df_db = pd.read_csv( experiment_home + spr.args.database+ spr.args.lr_db,sep='\t', usecols=['lr_pair'])

states = [10,7,4,2,22,24,18,1]
top_lr = 200
keep_db= True
_network.interaction_statewise_lr_network(spr,states,top_lr,keep_db,df_db)

top_lr = 25
keep_db=False
_network.interaction_statewise_lr_network(spr,states,top_lr,keep_db,df_db)

##################################################################
# 6 look at neighbours of cancer cells gene expression
'''
-for only cancer cells select its neighbour cell type and interaction topic
-get raw data for these cancer cells and nbr cells
- normalize rowwise such that each ligand expression is proportion of all ligand expression
- get top 25 ligand and receptor genes for interaction topic
- filter those top genes and take mean, add 1e-5, and apply sqrt for heatmap
'''


from analysis import _topics
import matplotlib.pylab as plt
import colorcet as cc
import seaborn as sns
plt.rcParams['figure.figsize'] = [15, 4]
plt.rcParams['figure.autolayout'] = True

df=pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics_cancer_cells.csv.gz',compression='gzip')

cslist = [ 
('Epithelial_Normal',10),
('Epithelial_Normal',4),
('B_cells',24),
('T_cells',22),
('B_cells',22),
('T_cells',24),
('B_cell_n',22),
('B_cell_n',24),
('T_cell_n',22),
('T_cell_n',24),
('PVL',18),
('Myeloid',24),
('Endothelial',2)
]
_topics.plt_cn(spr,df,cslist)


##################################################################
# 7 interaction topic gse analysis
##################################################################
'''
take ligand and receptor weight and average them
'''
df_db = pd.read_csv( experiment_home + spr.args.database+ spr.args.lr_db,sep='\t', usecols=['lr_pair'])
df = _gsea.gse_interactiontopic_v2(spr,df_db)
df.to_csv(spr.interaction_topic.model_id+'_it_gsea.csv.gz',index=False,compression='gzip')



from analysis import _gsea
_gsea.gsea(spr)


