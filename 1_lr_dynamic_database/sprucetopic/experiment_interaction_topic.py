import os 
import sys
import datetime
from util._io import read_config
from collections import namedtuple
import logging
import pandas as pd
import spruce
from analysis import _network
import torch

mode= sys.argv[1]
now = datetime.datetime.now()
# args_home = '/home/BCCRC.CA/ssubedi/projects/spruce_topic/'
args_home = '/home/sishirsubedi/projects/spruce_topic/'

# os.chdir(args_home)
os.environ['args_home'] = args_home

params = read_config(args_home+'config/bcmix.yaml')
args = namedtuple('Struct',params.keys())(*params.values())


if mode =='train':
	model_info = args_home+args.output+args.interaction_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']
	id = now.strftime('%Y%m%d%H%M')
	model_id = model_info+'_'+id


	logging.basicConfig(filename=model_id+'.log',
					format='%(asctime)s %(levelname)-8s %(message)s',
					level=logging.INFO,
					datefmt='%Y-%m-%d %H:%M:%S')

	sp = spruce.Spruce()
	sp.data.raw_l_data = args_home+args.input+args.raw_l_data
	sp.data.raw_r_data = args_home+args.input+args.raw_r_data
	sp.data.cell_topic_h = args_home+args.output+args.cell_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']+'_cell_topic_h.tsv.gz'
	sp.data.neighbour = args_home+args.output+args.interaction_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']+'_nbr.pkl'

	sp.model_id = model_id
	
	batch_size = args.interaction_topic['train']['batch_size']
	epochs = args.interaction_topic['train']['epochs']
	layers1 = args.interaction_topic['train']['layers1']
	layers2 = args.interaction_topic['train']['layers2']
	latent_dims = args.interaction_topic['train']['latent_dims']
	input_dims1 = args.interaction_topic['train']['input_dims1']
	input_dims2 = args.interaction_topic['train']['input_dims2']

	f_loss = sp.model_id + '_loss.txt'
	device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
	
	model = sp.run_interaction_topic_ddb(batch_size,epochs,layers1,layers2,latent_dims,input_dims1,input_dims2,device,f_loss)

	torch.save(model.state_dict(), sp.model_id + '_interaction_topic_ddb.torch')

elif mode=='eval':
	sp = spruce.Spruce()
	model_info = args_home+args.output+args.interaction_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']
	model_id = model_info+args.interaction_topic['model_info']
	sp.model_id = model_id

	sp.interaction_topic.model = torch.load(sp.model_id + '_interaction_topic_ddb.torch')

	layers1 = args.interaction_topic['train']['layers1']
	layers2 = args.interaction_topic['train']['layers2']
	latent_dims = args.interaction_topic['train']['latent_dims']
	input_dims1 = args.interaction_topic['train']['input_dims1']
	input_dims2 = args.interaction_topic['train']['input_dims2']
	device = 'cpu'

	# df_alpha,beta,df_alpha_bias,beta_bias = sp.eval_interaction_topic_ddb(input_dims1,input_dims2,latent_dims,layers1,layers2)
	# df_alpha.to_csv(sp.model_id+'_ietm_beta1.tsv.gz',sep='\t',index=False,compression='gzip')
	# df_alpha_bias.to_csv(sp.model_id+'_ietm_beta2.tsv.gz',sep='\t',index=False,compression='gzip')

	# for i in range(25):
	# 	df_beta = pd.DataFrame(beta[i])
	# 	df_beta.to_csv(sp.model_id+'topic_'+str(i)+'_ietm_beta.tsv.gz',sep='\t',index=False,compression='gzip')

	# for i in range(25):
	# 	df_beta_bias = pd.DataFrame(beta_bias[i])
	# 	df_beta_bias.to_csv(sp.model_id+'topic_'+str(i)+'_ietm_beta_bias.tsv.gz',sep='\t',index=False,compression='gzip')

	df_beta,alpha,df_beta_bias,alpha_bias = sp.eval_interaction_topic_ddb(input_dims1,input_dims2,latent_dims,layers1,layers2)
	df_beta.to_csv(sp.model_id+'_ietm_beta.tsv.gz',sep='\t',index=False,compression='gzip')
	df_beta_bias.to_csv(sp.model_id+'_ietm_beta_bias.tsv.gz',sep='\t',index=False,compression='gzip')

	for i in range(25):
		df_alpha = pd.DataFrame(alpha[i])
		df_alpha.to_csv(sp.model_id+'topic_'+str(i)+'_ietm_alpha.tsv.gz',sep='\t',index=False,compression='gzip')

	for i in range(25):
		df_alpha_bias = pd.DataFrame(alpha_bias[i])
		df_alpha_bias.to_csv(sp.model_id+'topic_'+str(i)+'_ietm_alpha_bias.tsv.gz',sep='\t',index=False,compression='gzip')

	
elif mode=='results':
	sp = spruce.Spruce()
	model_info = args_home+args.output+args.interaction_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']
	model_id = model_info+args.interaction_topic['model_info']
	sp.model_id = model_id

	# sp.data.raw_l_data_genes = [x[0] for x in pd.read_pickle(args_home+args.input+args.raw_l_data_genes).values]
	# sp.data.raw_r_data_genes = [x[0] for x in pd.read_pickle(args_home+args.input+args.raw_r_data_genes).values]
	# _network.lr_network(sp)

	import seaborn as sns
	import matplotlib.pylab as plt

	# for i in range(25):	
	# 	df = pd.read_csv(sp.model_id+'topic_'+str(i)+'_ietm_alpha.tsv.gz',sep='\t',compression='gzip')
	# 	sns.heatmap(df, cmap="PuRd")
	# 	plt.savefig(sp.model_id+'topic_'+str(i)+'_heatmap.png')
	# 	plt.close()

	df = pd.read_pickle(args_home+args.input+args.raw_lr_data)
	sns.heatmap(df, cmap="PuRd")
	plt.savefig('cell_talk_db_heatmap.png')
	plt.close()

