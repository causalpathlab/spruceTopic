import os 
import datetime
from _utils.io import read_config
from collections import namedtuple
import logging
import pandas as pd

now = datetime.datetime.now()
args_home = '/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/'

os.chdir(args_home)
os.environ['args_home'] = args_home

params = read_config(args_home+'/config/pbmc.yaml')
args = namedtuple('Struct',params.keys())(*params.values())
model_file = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+now.strftime('%Y%m%d%H%M')

print(model_file)
# logging.basicConfig(filename=model_file+'.log',
# 						format='%(asctime)s %(levelname)-8s %(message)s',
# 						level=logging.INFO,
# 						datefmt='%Y-%m-%d %H:%M:%S')

##
from plot import plt_umap 
plt_umap.plot_umap_from_model(args)

##
from evals import res_pbmcs
res_pbmcs.pbmc_sample_cells_with_latent(args)

##
from evals import results
results.topic_top_genes(args)