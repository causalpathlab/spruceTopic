import os
import sys
import datetime
from util._io import read_config
from collections import namedtuple
from pathlib import Path
from model import _cell_topic
import logging
import pandas as pd
import spruce
import torch

def run_model(experiment_home,args,mode):

	now = datetime.datetime.now()

	if mode =='train':
		model_info = experiment_home+args.output+args.cell_topic['out']+args.cell_topic['model_info']
		id = now.strftime('%Y%m%d')
		model_id = model_info+'_'+id


		logging.basicConfig(filename=model_id+'.log',
						format='%(asctime)s %(levelname)-8s %(message)s',
						level=logging.INFO,
						datefmt='%Y-%m-%d %H:%M:%S')

		sp = spruce.Spruce()
		sp.data.sparse_data = experiment_home + args.input + args.sparse_data 
		sp.data.sparse_data_ids = experiment_home + args.input + args.sparse_data_ids
		sp.model_id = model_id

		batch_size = args.cell_topic['train']['batch_size']
		l_rate = args.cell_topic['train']['l_rate']
		epochs = args.cell_topic['train']['epochs']
		layers = args.cell_topic['train']['layers']
		latent_dims = args.cell_topic['train']['latent_dims']
		device = args.cell_topic['train']['device']

		loss_values = sp.run_cell_topic(batch_size,l_rate,epochs,layers,latent_dims,device)
		torch.save(sp.cell_topic.model.state_dict(), sp.model_id +  args.fn_cell_topic_model)
		dflv = pd.DataFrame(loss_values[0])
		dflv.to_csv(sp.model_id + args.fn_cell_topic_loss,index=False,compression='gzip')
		dflv = pd.DataFrame(loss_values[1])
		dflv.to_csv(sp.model_id + args.fn_cell_topic_loss,index=False,compression='gzip')

	elif mode=='eval':
		sp = spruce.Spruce()
		sp.data.sparse_data = experiment_home + args.input + args.sparse_data 
		sp.data.sparse_data_ids = experiment_home + args.input + args.sparse_data_ids
		sp.model_id = experiment_home+args.output+args.cell_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']

		sp.cell_topic.model = torch.load(sp.model_id + args.fn_cell_topic_model)
		batch_size = args.cell_topic['eval']['batch_size']
		layers = args.cell_topic['train']['layers']
		latent_dims = args.cell_topic['train']['latent_dims']
		device = 'cpu'
		df_z,df_h,df_beta,df_beta_var = sp.eval_cell_topic(batch_size,layers,latent_dims,device)
		df_z.to_csv(sp.model_id+ args.fn_cell_topic_z,index=False,compression='gzip')
		df_h.to_csv(sp.model_id+ args.fn_cell_topic_h,index=False,compression='gzip')
		df_beta.to_csv(sp.model_id+ args.fn_cell_topic_beta_mean,index=False,compression='gzip')
		df_beta_var.to_csv(sp.model_id + args.fn_cell_topic_beta_var,index=False,compression='gzip')

def get_model(experiment_home,args):

	sp = spruce.Spruce()
	sp.model_id = experiment_home+args.cell_topic['out']+args.cell_topic['model_info']+args.cell_topic['model_id']

	input_dims = args.cell_topic['eval']['input_dims']
	layers = args.cell_topic['train']['layers']
	latent_dims = args.cell_topic['train']['latent_dims']
	device = 'cpu'

	model = _cell_topic.ETM(input_dims,latent_dims,layers).to(device)
	model.load_state_dict(torch.load(sp.model_id + args.fn_cell_topic_model,map_location=torch.device('cpu')))
	model.eval()
	sp.cell_topic.model = model

	sp.cell_topic.z = pd.read_csv(sp.model_id + args.fn_cell_topic_z,compression='gzip')
	sp.cell_topic.h = pd.read_csv(sp.model_id + args.fn_cell_topic_h,compression='gzip')
	
	sp.data.raw_data_genes = pd.read_pickle(experiment_home+args.input+args.raw_data_genes)[0].values

	sp.cell_topic.beta_mean = pd.read_csv(sp.model_id+ args.fn_cell_topic_beta_mean,compression='gzip')
	sp.cell_topic.beta_var = pd.read_csv(sp.model_id+ args.fn_cell_topic_beta_var,compression='gzip')

	sp.cell_topic.beta_mean.columns = sp.data.raw_data_genes
	sp.cell_topic.beta_var.columns = sp.data.raw_data_genes

	
	return sp

def get_experiment(experiment):

	server = Path.home().as_posix()
	experiment_home = server+experiment
	experiment_config = read_config(experiment_home+'config.yaml')
	args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())

	return get_model(experiment_home,args)


# def temp(sp):
# 	print('processing...',sp.model_id)
# 	# _umap_viz.plot_umap_with_annotation_mix(sp)
# 	# _topics.topic_top_genes(sp,5)
# 	# _topics.topic_top_genes(sp,10)
# 	# _topics.topic_top_genes(sp,25)
# 	# _topics.sample_cells_with_latent(sp)
# 	# _umap_viz.plot_umap(sp)


# 	hallmark_db = args.home + args.database +  args.db_hallmark
# 	df_hmark = pd.read_csv(hallmark_db,compression='gzip')
# 	df_hmark = df_hmark[['gs_name','gene_symbol']]
# 	df_htest = _gsea.hypergeom_test(sp,df_hmark)
# 	df_htest = df_htest.sort_values('pval')
# 	df_htest.to_csv(sp.model_id+'_hallmark_hypergeom_test.tsv.gz',compression='gzip',index=False,sep='\t')	

# 	hallmark_db = args.home + args.database +  args.db_cmark
# 	df_hmark = pd.read_excel(hallmark_db)
# 	df_hmark = df_hmark[['Cell type','Gene']]
# 	df_hmark.columns = ['gs_name','gene_symbol']
# 	df_htest = _gsea.hypergeom_test(sp,df_hmark)
# 	df_htest = df_htest.sort_values('pval')
# 	df_htest.to_csv(sp.model_id+'_cmark_hypergeom_test.tsv.gz',compression='gzip',index=False,sep='\t')	
	
	




if __name__ == "__main__":
	

	experiment = sys.argv[1]
	mode = sys.argv[2]

	server = Path.home().as_posix()
	experiment_home = server+experiment
	experiment_config = read_config(experiment_home+'config.yaml')
	args = namedtuple('Struct',experiment_config.keys())(*experiment_config.values())
	
	run_model(experiment_home,args,mode)