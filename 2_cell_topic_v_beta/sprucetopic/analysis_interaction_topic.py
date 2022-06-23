
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
    # analysis of latent h
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

# df_its.to_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz',index=False,compression='gzip')
# df_its = pd.read_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz',compression='gzip')

df_hmax = pd.DataFrame(pd.Series(df_its.iloc[:,1:].values.flatten()).value_counts()).reset_index().rename(columns={'index':'topic',0:'argmax_count'})
df_hmax = df_hmax.sort_values('argmax_count',ascending=False)
p = sns.barplot(x='topic',y='argmax_count',data=df_hmax,color='blue')
p.set_xlabel("Topic",fontsize=20)
p.set_ylabel("Count(argmax)",fontsize=20)
plt.savefig(spr.interaction_topic.model_id+'_it_h_argmax.png');plt.close()


####### generate topics summary 
'''
Take umap coordinates, topic assignment, cluster assignment, cell type, and cluster majority cell type from
cell topic analysis and merge with state assignment from interaction topic.
'''
spr.interaction_topic.neighbour_h =  pd.read_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz')
df_kmeans = pd.read_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_kmeans.csv.gz')
df_all_topics = _topics.get_topics(spr,df_kmeans)
df_all_topics.to_csv(spr.interaction_topic.model_id+'_it_model_all_topics.csv.gz',index=False,compression='gzip')

'''
Generate a model summary file where each cluster is assigned-
1. cell topic based on argmax
2. interaction state based on highest state assigned among all neighbours
3. total number of cells in the cluster
'''
dfsummary = _topics.topics_summary(spr,df_kmeans)	
### filter topics summary for summary plot
clust_to_cell_type = dict(zip( df_kmeans['cluster'],df_kmeans['cluster_celltype']))
dfsummary['celltype'] = [clust_to_cell_type[x] for x in dfsummary['cluster']]
dfsummary = dfsummary.sort_values('state')
min_cells_per_state = 100
dfsummary = dfsummary[dfsummary['cell']>min_cells_per_state]
dfsummary.to_csv(spr.interaction_topic.model_id+'_it_model_summary.csv.gz',index=False,compression='gzip')


##################################################################
    # network graph
##################################################################

from analysis import _network
df = pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics.csv.gz',compression='gzip')
_network.cell_interaction_network(spr,df)


##################################################################
    # top genes
##################################################################


from analysis import _topics

#### get top genes based on beta weight matrix - l/r separately

top_n = 10
df_top_genes = _topics.topic_top_lr_genes(spr,top_n)
df_top_genes.to_csv(spr.interaction_topic.model_id+'_it_beta_weight_top_'+str(top_n)+'_genes.csv.gz',index=False,compression='gzip')

#### get top genes based on beta weight matrix - l/r separately

df_db = pd.read_csv( experiment_home + spr.args.database+ spr.args.lr_db,sep='\t', usecols=['lr_pair'])
top_n=25
df_lr_topic = _topics.topic_top_lr_pair_genes(spr,df_db,top_n)
df_lr_topic.to_csv(spr.interaction_topic.model_id+'_it_beta_top_'+str(top_n)+'_lrpair.csv.gz',compression='gzip')

##################################################################
    # cancer centric view
##################################################################

from analysis import _topics
df = pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics.csv.gz',compression='gzip')
df_its = pd.read_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz',compression='gzip')
df_nbr = spr.cell_topic.neighbour

res = _topics.get_cell_neighbours_states(df_nbr,df_its,df)
df_res = pd.DataFrame(res,columns=['cancer_cell','nbr','state'])
df_res = df_res.explode(['nbr','state'])
df_res = pd.merge(df_res,df[['cell','cluster','cluster_celltype']],left_on='nbr',right_on='cell',how='left')
df_res.to_csv(spr.interaction_topic.model_id+'_it_model_all_topics_cancer_cells.csv.gz',index=False,compression='gzip')


##### look at non cancer neighbours of cancer cells and check state distribution 
df=pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics_cancer_cells.csv.gz',compression='gzip')

df = df[df['cluster_celltype'] != 'Cancer Epithelial']
df_grp = df.groupby(['cluster_celltype','state'])['state'].size().rename('count').reset_index()
##normalize
celltype_sum = dict(df_grp.groupby('cluster_celltype')['count'].sum())
df_grp['ncount'] = [x/celltype_sum[y] for x,y in zip(df_grp['count'],df_grp['cluster_celltype'])]
df_grp.to_csv(spr.interaction_topic.model_id+'_it_cancercells_interactions.csv.gz',index=False,compression='gzip')

### look at cancer neighbours of cancer cells


from analysis import _topics
import matplotlib.pylab as plt
import colorcet as cc
import seaborn as sns
plt.rcParams['figure.figsize'] = [15, 4]
plt.rcParams['figure.autolayout'] = True

df=pd.read_csv(spr.interaction_topic.model_id+'_it_model_all_topics_cancer_cells.csv.gz',compression='gzip')

cslist = [ 
('Normal Epithelial',10),
('Normal Epithelial',4),
('B-cells',22),
('T-cells',22),
('PVL',18),
('Myeloid',24),
('Endothelial',2)
]
_topics.plt_cn(spr,df,cslist)