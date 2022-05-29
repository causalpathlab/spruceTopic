import os 
import sys
import datetime
from util._io import read_config
from collections import namedtuple
import logging
import pandas as pd
import spruce
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
	model_info = args_home+args.output+args.cell_topic['out']+args.cell_topic['model_info']
	id = now.strftime('%Y%m%d%H%M')
	model_id = model_info+'_'+id


	logging.basicConfig(filename=model_id+'.log',
					format='%(asctime)s %(levelname)-8s %(message)s',
					level=logging.INFO,
					datefmt='%Y-%m-%d %H:%M:%S')

	spruce = spruce.Spruce()
	spruce.data.sparse_data = args_home + args.input + args.sparse_data 
	spruce.data.sparse_data_ids = args_home + args.input + args.sparse_data_ids
	spruce.model_id = model_id

	batch_size = args.cell_topic['train']['batch_size']
	l_rate = args.cell_topic['train']['l_rate']
	epochs = args.cell_topic['train']['epochs']
	layers = args.cell_topic['train']['layers']
	latent_dims = args.cell_topic['train']['latent_dims']
	device = args.cell_topic['train']['device']

	loss_values = spruce.run_cell_topic_model(batch_size,l_rate,epochs,layers,latent_dims,device)
	torch.save(spruce.cell_topic.model.state_dict(), spruce.model_id + '_cell_topic.torch')
	dflv = pd.DataFrame(loss_values[0])
	dflv.to_csv(spruce.model_id + '_cell_topic_loss.txt.gz',index=False,compression='gzip')
	dflv = pd.DataFrame(loss_values[1])
	dflv.to_csv(spruce.model_id + '_cell_topic_loss2.txt.gz',index=False,compression='gzip')

elif mode=='eval':
	spruce = spruce.Spruce()
	spruce.data.sparse_data = args_home + args.input + args.sparse_data 
	spruce.data.sparse_data_ids = args_home + args.input + args.sparse_data_ids
	model_info = args_home+args.output+args.cell_topic['out']+args.cell_topic['model_info']
	id = '202205290330'
	model_id = model_info+'_'+id
	spruce.model_id = model_id

	spruce.cell_topic.model = torch.load(spruce.model_id + '_cell_topic.torch')
	batch_size = args.cell_topic['eval']['batch_size']
	layers = args.cell_topic['train']['layers']
	latent_dims = args.cell_topic['train']['latent_dims']
	device = 'cpu'
	df_z,df_h,df_beta1 = spruce.eval_cell_topic_model(batch_size,layers,latent_dims,device)
	df_z.to_csv(spruce.model_id+'_cell_topic_z.tsv.gz',sep='\t',index=False,compression='gzip')
	df_h.to_csv(spruce.model_id+'_cell_topic_h.tsv.gz',sep='\t',index=False,compression='gzip')
	df_beta1.to_csv(spruce.model_id+'_cell_topic_beta.tsv.gz',sep='\t',index=False,compression='gzip')