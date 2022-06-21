
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


df_its = spr.interaction_topic_states()
# df_its.to_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz',index=False,compression='gzip')

df_hmax = pd.DataFrame(pd.Series(df_its.iloc[:,1:].values.flatten()).value_counts()).reset_index().rename(columns={'index':'topic',0:'argmax_count'})
df_hmax = df_hmax.sort_values('argmax_count',ascending=False)
p = sns.barplot(x='topic',y='argmax_count',data=df_hmax,color='blue')
p.set_xlabel("Topic",fontsize=20)
p.set_ylabel("Count(argmax)",fontsize=20)
plt.savefig(spr.interaction_topic.model_id+'_it_h_argmax.png');plt.close()


####### generate topics summary 

spr.interaction_topic.neighbour_h =  pd.read_csv(spr.interaction_topic.model_id+'_it_h_argmax.csv.gz')
df_kmeans = pd.read_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_kmeans.csv.gz')
df_all_topics = _topics.get_topics(spr,df_kmeans)
df_all_topics.to_csv(spr.interaction_topic.model_id+'_it_model_all_topics.csv.gz',index=False,compression='gzip')


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

top_n = 5
df_top_genes = _topics.topic_top_lr_genes(spr,top_n)
df_top_genes.to_csv(spr.interaction_topic.model_id+'_it_beta_weight_top_'+str(top_n)+'_genes.csv.gz',index=False,compression='gzip')

#### get top genes based on beta weight matrix - l/r separately

df_db = pd.read_csv( experiment_home + spr.args.database+ spr.args.lr_db,sep='\t', usecols=['lr_pair'])
top_n=25
df_lr_topic = _topics.topic_top_lr_pair_genes(spr,df_db,top_n)
df_lr_topic.to_csv(spr.interaction_topic.model_id+'_it_beta_top_'+str(top_n)+'_lrpair.csv.gz',compression='gzip')

