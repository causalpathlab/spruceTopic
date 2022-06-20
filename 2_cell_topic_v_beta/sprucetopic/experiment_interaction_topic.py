import os 
import sys
import datetime
from util._io import read_config
from collections import namedtuple
from pathlib import Path
import logging
import pandas as pd
import spruce
from analysis import _topics
from model import _interaction_topic,_neighbour
import torch


def run_model(experiment_home,args,mode):

	if mode == 'nbr':
		
		# id = datetime.datetime.now().strftime('%Y%m%d')

		model_id = experiment_home+ args.interaction_topic['out']+args.cell_topic['model_id']

		logging.basicConfig(filename=model_id+'_generate_nbr.log',
						format='%(asctime)s %(levelname)-8s %(message)s',
						level=logging.INFO,
						datefmt='%Y-%m-%d %H:%M:%S')

		sp = spruce.Spruce()
		sp.args = args
		sp.model_id = model_id

		sp.data.raw_l_data_genes = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_l_data_genes)[0].values
		sp.data.raw_r_data_genes = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_r_data_genes)[0].values
		sp.data.raw_l_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_l_data)
		sp.data.raw_r_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_r_data)
		sp.data.raw_lr_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_lr_data)

		sp.cell_topic.h = pd.read_csv(experiment_home+args.cell_topic['out']+args.cell_topic['model_id']+'_ct_h.csv.gz',compression='gzip')

		df_nbr = _neighbour.generate_neighbours(sp)
		df_nbr = df_nbr.rename(columns={'Unnamed: 0':'cell'})
		df_nbr.to_csv(sp.model_id+'_nbr_cellids.csv.gz',compression='gzip')



	elif mode =='train':

		# id = datetime.datetime.now().strftime('%Y%m%d%H')

		# model_dir = experiment_home+args.interaction_topic['out']+id
		
		# args.cell_topic['model_info']+args.cell_topic['model_id']
		# model_id = 'model'

		model_id = experiment_home+ args.interaction_topic['out']+args.cell_topic['model_id']

		logging.basicConfig(filename=model_id+'_it_train.log',
						format='%(asctime)s %(levelname)-8s %(message)s',
						level=logging.INFO,
						datefmt='%Y-%m-%d %H:%M:%S')

		sp = spruce.Spruce()
		sp.data.raw_l_data = experiment_home+args.data+args.sample_id+args.raw_l_data
		sp.data.raw_r_data = experiment_home+args.data+args.sample_id+args.raw_r_data
		sp.data.raw_lr_data = experiment_home+args.data+args.sample_id+args.raw_lr_data
		sp.cell_topic.h = experiment_home+args.cell_topic['out']+args.cell_topic['model_id']+'_ct_h.csv.gz'
		sp.data.neighbour = experiment_home+args.cell_topic['out']+args.cell_topic['model_id']+'_nbr_cellids.csv.gz'

		sp.model_id = model_id
		
		batch_size = args.interaction_topic['train']['batch_size']
		epochs = args.interaction_topic['train']['epochs']
		layers1 = args.interaction_topic['train']['layers1']
		layers2 = args.interaction_topic['train']['layers2']
		latent_dims = args.interaction_topic['train']['latent_dims']
		input_dims1 = args.interaction_topic['train']['input_dims1']
		input_dims2 = args.interaction_topic['train']['input_dims2']

		f_loss = sp.model_id + '_it_model_loss.txt'
		device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
		
		model = sp.run_interaction_topic(batch_size,epochs,layers1,layers2,latent_dims,input_dims1,input_dims2,device,f_loss)

		torch.save(model.state_dict(), sp.model_id + '_it_model.torch')

		df = pd.read_csv(f_loss,header=None,sep='\t')
		df.groupby(df.index//4500).mean().to_csv(sp.model_id+'_it_lossm.txt.gz',compression='gzip',index=False,header=None)	



	elif mode=='eval':
		sp = spruce.Spruce()
		sp.model_id = experiment_home+args.interaction_topic['out']+args.interaction_topic['model_id']

		sp.interaction_topic.model = torch.load(sp.model_id + '_it_model.torch')
		batch_size = args.interaction_topic['train']['batch_size']

		layers1 = args.interaction_topic['train']['layers1']
		layers2 = args.interaction_topic['train']['layers2']
		latent_dims = args.interaction_topic['train']['latent_dims']
		input_dims1 = args.interaction_topic['train']['input_dims1']
		input_dims2 = args.interaction_topic['train']['input_dims2']
		device = 'cpu'

		df_betal,df_betar,df_betal_bias,df_betar_bias = sp.eval_interaction_topic(batch_size,input_dims1,input_dims2,latent_dims,layers1,layers2)

		df_betal.to_csv(sp.model_id+'_it_beta_l.tsv.gz',sep='\t',index=False,compression='gzip')
		df_betar.to_csv(sp.model_id+'_it_beta_r.tsv.gz',sep='\t',index=False,compression='gzip')
		df_betal_bias.to_csv(sp.model_id+'_it_beta_l_bias.tsv.gz',sep='\t',index=False,compression='gzip')
		df_betar_bias.to_csv(sp.model_id+'_it_beta_r_bias.tsv.gz',sep='\t',index=False,compression='gzip')



def get_model(experiment_home,args):

	sp = spruce.Spruce()
	sp.args = args
	sp.model_id = experiment_home+ args.interaction_topic['out']+args.interaction_topic['model_id']

	sp.data.raw_l_data_genes = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_l_data_genes)[0].values
	sp.data.raw_r_data_genes = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_r_data_genes)[0].values
	sp.data.raw_l_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_l_data)
	sp.data.raw_r_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_r_data)
	sp.data.raw_lr_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_lr_data)

	sp.cell_topic.h = pd.read_csv(experiment_home+args.cell_topic['out']+args.cell_topic['model_id']+'_ct_h.csv.gz')
	sp.data.neighbour = pd.read_csv(experiment_home+args.cell_topic['out']+args.cell_topic['model_id']+'_nbr_cellids.csv.gz')

	batch_size = sp.args.interaction_topic['train']['batch_size']
	layers1 = sp.args.interaction_topic['train']['layers1']
	layers2 = sp.args.interaction_topic['train']['layers2']
	latent_dims = sp.args.interaction_topic['train']['latent_dims']
	input_dims1 = sp.args.interaction_topic['train']['input_dims1']
	input_dims2 = sp.args.interaction_topic['train']['input_dims2']

	model = _interaction_topic.LitETM(batch_size,input_dims1,input_dims2,latent_dims,layers1,layers2,'_txt_')
	model.load_state_dict(torch.load(sp.model_id+'_it_model.torch'))
	model.eval()
	sp.interaction_topic.model = model

	sp.interaction_topic.beta_l = pd.read_csv(sp.model_id+'_it_beta_l.tsv.gz',sep='\t',compression='gzip')
	sp.interaction_topic.beta_r = pd.read_csv(sp.model_id+'_it_beta_r.tsv.gz',sep='\t',compression='gzip')

	sp.interaction_topic.beta_l.columns = sp.data.raw_r_data_genes
	sp.interaction_topic.beta_r.columns = sp.data.raw_l_data_genes

	return sp


def get_experiment_model(experiment):

	server = Path.home().as_posix()
	experiment_home = server+experiment
	experiment_config = read_config(experiment_home+'config.yaml')
	args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())
	return get_model(experiment_home,args)

if __name__ == "__main__":
	

	experiment = sys.argv[1]
	mode = sys.argv[2]

	server = Path.home().as_posix()
	experiment_home = server+experiment
	experiment_config = read_config(experiment_home+'config.yaml')
	args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())
	
	run_model(experiment_home,args,mode)