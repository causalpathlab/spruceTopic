

import experiment_interaction_topic 
import pandas as pd


experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'
spr = experiment_interaction_topic.get_experiment_model(experiment_home)
print(spr.model_id)

##################################################################
    # top genes
##################################################################


from analysis import _topics

top_n = 5
df_top_genes = _topics.topic_top_lr_genes(spr,top_n)
df_top_genes.to_csv(spr.model_id+'_it_beta_weight_top_'+str(top_n)+'_genes.csv.gz',index=False,compression='gzip')

from pathlib import Path
server = Path.home().as_posix()
experiment_home = server+experiment_home
df_db = pd.read_csv( experiment_home + spr.args.database+ spr.args.lr_db,sep='\t', usecols=['lr_pair'])
top_n=25
df_lr_topic = _topics.topic_top_lr_pair_genes(spr,df_db,top_n)
df_lr_topic.to_csv(spr.model_id+'_it_beta_top_'+str(top_n)+'_lrpair.csv.gz',compression='gzip')

