import os 
import sys
import datetime
from util._io import read_config
from collections import namedtuple
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
		sp.interaction_topic.model_id = model_id

		sp.data.raw_l_data_genes = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_l_data_genes)[0].values
		sp.data.raw_r_data_genes = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_r_data_genes)[0].values
		sp.data.raw_l_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_l_data)
		sp.data.raw_r_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_r_data)
		sp.data.raw_lr_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_lr_data)

		sp.cell_topic.h = pd.read_csv(experiment_home+args.cell_topic['out']+args.cell_topic['model_id']+'_ct_h.csv.gz',compression='gzip')

		df_nbr = _neighbour.generate_neighbours(sp)
		df_nbr = df_nbr.rename(columns={'Unnamed: 0':'cell'})
		df_nbr.to_csv(sp.interaction_topic.model_id+'_nbr_cellids.csv.gz',compression='gzip')



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
		sp.cell_topic.neighbour = experiment_home+args.cell_topic['out']+args.cell_topic['model_id']+'_nbr_cellids.csv.gz'

		sp.interaction_topic.model_id = model_id
		
		batch_size = args.interaction_topic['train']['batch_size']
		epochs = args.interaction_topic['train']['epochs']
		layers1 = args.interaction_topic['train']['layers1']
		layers2 = args.interaction_topic['train']['layers2']
		latent_dims = args.interaction_topic['train']['latent_dims']
		input_dims1 = args.interaction_topic['train']['input_dims1']
		input_dims2 = args.interaction_topic['train']['input_dims2']

		f_loss = sp.interaction_topic.model_id + '_it_model_loss.txt'
		device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
		
		model = sp.run_interaction_topic(batch_size,epochs,layers1,layers2,latent_dims,input_dims1,input_dims2,device,f_loss)

		torch.save(model.state_dict(), sp.interaction_topic.model_id + '_it_model.torch')

		df = pd.read_csv(f_loss,header=None,sep='\t')
		df.groupby(df.index//4500).mean().to_csv(sp.interaction_topic.model_id+'_it_model_lossm.txt.gz',compression='gzip',index=False,header=None)	



	elif mode=='eval':
		sp = spruce.Spruce()
		sp.interaction_topic.model_id = experiment_home+args.interaction_topic['out']+args.interaction_topic['model_id']

		sp.interaction_topic.model = torch.load(sp.interaction_topic.model_id + '_it_model.torch')
		batch_size = args.interaction_topic['train']['batch_size']

		layers1 = args.interaction_topic['train']['layers1']
		layers2 = args.interaction_topic['train']['layers2']
		latent_dims = args.interaction_topic['train']['latent_dims']
		input_dims1 = args.interaction_topic['train']['input_dims1']
		input_dims2 = args.interaction_topic['train']['input_dims2']
		device = 'cpu'

		df_betalm,df_betalv,df_betarm,df_betarv,df_betal_bias,df_betar_bias = sp.eval_interaction_topic(batch_size,input_dims1,input_dims2,latent_dims,layers1,layers2)

		df_betalm.to_csv(sp.interaction_topic.model_id+'_it_beta_lm.csv.gz',index=False,compression='gzip')
		df_betalv.to_csv(sp.interaction_topic.model_id+'_it_beta_lv.csv.gz',index=False,compression='gzip')
		df_betarm.to_csv(sp.interaction_topic.model_id+'_it_beta_rm.csv.gz',index=False,compression='gzip')
		df_betarv.to_csv(sp.interaction_topic.model_id+'_it_beta_rv.csv.gz',index=False,compression='gzip')
		df_betal_bias.to_csv(sp.interaction_topic.model_id+'_it_beta_l_bias.csv.gz',index=False,compression='gzip')
		df_betar_bias.to_csv(sp.interaction_topic.model_id+'_it_beta_r_bias.csv.gz',index=False,compression='gzip')

def get_model(experiment_home,args):

	sp = spruce.Spruce()
	sp.args = args
	sp.cell_topic.model_id = experiment_home+ args.cell_topic['out']+args.cell_topic['model_id']
	sp.cell_topic.id = experiment_home+ args.cell_topic['out']

	sp.interaction_topic.model_id = experiment_home+ args.interaction_topic['out']+args.interaction_topic['model_id']
	sp.interaction_topic.id = experiment_home+ args.interaction_topic['out']

	sp.data.raw_l_data_genes = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_l_data_genes)[0].values
	sp.data.raw_r_data_genes = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_r_data_genes)[0].values
	sp.data.raw_l_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_l_data)
	sp.data.raw_r_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_r_data)
	sp.data.raw_lr_data = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_lr_data)

	sp.cell_topic.h = pd.read_csv(experiment_home+args.cell_topic['out']+args.cell_topic['model_id']+'_ct_h.csv.gz')
	sp.cell_topic.beta_mean = pd.read_csv(experiment_home+args.cell_topic['out']+args.cell_topic['model_id']+'_ct_beta_mean.csv.gz')
	sp.data.raw_data_genes = pd.read_pickle(experiment_home+args.data+args.sample_id+args.raw_data_genes)[0].values
	sp.cell_topic.beta_mean.columns = sp.data.raw_data_genes


	sp.cell_topic.neighbour = pd.read_csv(experiment_home+args.cell_topic['out']+args.cell_topic['model_id']+'_nbr_cellids.csv.gz')

	batch_size = sp.args.interaction_topic['train']['batch_size']
	layers1 = sp.args.interaction_topic['train']['layers1']
	layers2 = sp.args.interaction_topic['train']['layers2']
	latent_dims = sp.args.interaction_topic['train']['latent_dims']
	input_dims1 = sp.args.interaction_topic['train']['input_dims1']
	input_dims2 = sp.args.interaction_topic['train']['input_dims2']

	model = _interaction_topic.LitETM(batch_size,input_dims1,input_dims2,latent_dims,layers1,layers2,'_txt_')
	model.load_state_dict(torch.load(sp.interaction_topic.model_id+'_it_model.torch'))
	model.eval()
	sp.interaction_topic.model = model

	sp.interaction_topic.beta_rm = pd.read_csv(sp.interaction_topic.model_id+'_it_beta_lm.csv.gz',compression='gzip')
	sp.interaction_topic.beta_rv = pd.read_csv(sp.interaction_topic.model_id+'_it_beta_lv.csv.gz',compression='gzip')
	sp.interaction_topic.beta_lm = pd.read_csv(sp.interaction_topic.model_id+'_it_beta_rm.csv.gz',compression='gzip')
	sp.interaction_topic.beta_lv = pd.read_csv(sp.interaction_topic.model_id+'_it_beta_rv.csv.gz',compression='gzip')

	sp.interaction_topic.beta_rm.columns = sp.data.raw_r_data_genes
	sp.interaction_topic.beta_rv.columns = sp.data.raw_r_data_genes
	sp.interaction_topic.beta_lm.columns = sp.data.raw_l_data_genes
	sp.interaction_topic.beta_lv.columns = sp.data.raw_l_data_genes

	return sp


def get_experiment_model(experiment_home):

	experiment_config = read_config(experiment_home+'config.yaml')
	args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())
	return get_model(experiment_home,args)

if __name__ == "__main__":
	

	experiment_home = sys.argv[1]
	mode = sys.argv[2]

	experiment_config = read_config(experiment_home+'config.yaml')
	args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())
	
	run_model(experiment_home,args,mode)