import os 
import sys
import datetime
from util._io import read_config
from collections import namedtuple
import logging
import pandas as pd
import spruce
from analysis import _topics
import torch

mode= sys.argv[1]
now = datetime.datetime.now()
# args_home = '/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/5_augmented_lr_multinomial_dirichlet_v_beta/'
args_home = '/home/sishirsubedi/projects/experiments/spruce_topic/5_augmented_lr_multinomial_dirichlet_v_beta/'

# os.chdir(args_home)
os.environ['args_home'] = args_home

params = read_config(args_home+'config/bcmix.yaml')
args = namedtuple('Struct',params.keys())(*params.values())


if mode =='train':
	model_info = args_home+args.output+args.interaction_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']
	id = now.strftime('%Y%m%d%H')
	model_id = model_info+'_'+id


	logging.basicConfig(filename=model_id+'.log',
					format='%(asctime)s %(levelname)-8s %(message)s',
					level=logging.INFO,
					datefmt='%Y-%m-%d %H:%M:%S')

	sp = spruce.Spruce()
	sp.data.raw_l_data = args.home+args.input+args.raw_l_data
	sp.data.raw_r_data = args.home+args.input+args.raw_r_data
	sp.data.raw_lr_data = args.home+args.input+args.raw_lr_data
	sp.cell_topic.h = args_home+args.output+args.cell_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']+'_cell_topic_h.tsv.gz'
	sp.data.neighbour = args_home+args.output+args.cell_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']+'_nbr.pkl'

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
	
	model = sp.run_interaction_topic(batch_size,epochs,layers1,layers2,latent_dims,input_dims1,input_dims2,device,f_loss)

	torch.save(model.state_dict(), sp.model_id + '_interaction_topic.torch')

elif mode=='eval':
	sp = spruce.Spruce()
	sp.model_id = args_home+args.output+args.interaction_topic['out']+args.interaction_topic['model_info']

	sp.interaction_topic.model = torch.load(sp.model_id + '_interaction_topic.torch')
	batch_size = args.interaction_topic['train']['batch_size']

	layers1 = args.interaction_topic['train']['layers1']
	layers2 = args.interaction_topic['train']['layers2']
	latent_dims = args.interaction_topic['train']['latent_dims']
	input_dims1 = args.interaction_topic['train']['input_dims1']
	input_dims2 = args.interaction_topic['train']['input_dims2']
	device = 'cpu'

	df_beta1,df_beta2,df_beta1_bias,df_beta2_bias = sp.eval_interaction_topic(batch_size,input_dims1,input_dims2,latent_dims,layers1,layers2)

	df_beta1.to_csv(sp.model_id+'_ietm_beta1.tsv.gz',sep='\t',index=False,compression='gzip')
	df_beta2.to_csv(sp.model_id+'_ietm_beta2.tsv.gz',sep='\t',index=False,compression='gzip')
	df_beta1_bias.to_csv(sp.model_id+'_ietm_beta1_bias.tsv.gz',sep='\t',index=False,compression='gzip')
	df_beta2_bias.to_csv(sp.model_id+'_ietm_beta2_bias.tsv.gz',sep='\t',index=False,compression='gzip')


elif mode=='results':
	sp = spruce.Spruce()
	model_info = args_home+args.output+args.interaction_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']
	id = '2022061216'
	model_id = model_info+'_'+id
	
	sp.model_id = model_id
	sp.interaction_topic.model = torch.load(sp.model_id + '_interaction_topic.torch')

	sp.data.raw_l_data = pd.read_pickle(args.home+args.input+args.raw_l_data)
	sp.data.raw_r_data = pd.read_pickle(args.home+args.input+args.raw_r_data)
	sp.data.raw_lr_data = pd.read_pickle(args.home+args.input+args.raw_lr_data)
	sp.cell_topic.h = pd.read_csv(args_home+args.output+args.cell_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']+'_cell_topic_h.tsv.gz',sep='\t',compression='gzip')
	sp.data.neighbour = pd.read_pickle(args_home+args.output+args.cell_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']+'_nbr.pkl')

	batch_size = args.interaction_topic['train']['batch_size']
	layers1 = args.interaction_topic['train']['layers1']
	layers2 = args.interaction_topic['train']['layers2']
	latent_dims = args.interaction_topic['train']['latent_dims']
	input_dims1 = args.interaction_topic['train']['input_dims1']
	input_dims2 = args.interaction_topic['train']['input_dims2']
	device = 'cpu'

	# topics_prob = sp.interactions_prob(input_dims1,input_dims2,latent_dims,layers1,layers2)
	# _topics.interaction_summary(sp,topics_prob)

	sp.interactions_state_summary(batch_size,input_dims1,input_dims2,latent_dims,layers1,layers2)


elif mode=='plots':

	sp = spruce.Spruce()
	model_info = args_home+args.output+args.interaction_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']
	id = '2022061216'
	model_id = model_info+'_'+id
	
	sp.model_id = model_id

	sp.data.raw_l_data_genes = pd.read_pickle(args.home+args.input+args.raw_l_data_genes)[0].values
	sp.data.raw_r_data_genes = pd.read_pickle(args.home+args.input+args.raw_r_data_genes)[0].values

	sp.interaction_topic.beta1 = pd.read_csv(sp.model_id+'_ietm_beta1.tsv.gz',sep='\t',compression='gzip')
	sp.interaction_topic.beta2 = pd.read_csv(sp.model_id+'_ietm_beta2.tsv.gz',sep='\t',compression='gzip')

	_topics.topic_top_lr_genes(sp)