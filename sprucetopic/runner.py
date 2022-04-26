import os 
import sys
import datetime
from _utils.io import read_config
from collections import namedtuple
import logging
import pandas as pd

mode= sys.argv[1]
now = datetime.datetime.now()
args_home = '/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/'
# args_home = '/home/sishirsubedi/projects/tumour_immune_interaction/'

# os.chdir(args_home)
os.environ['args_home'] = args_home

params = read_config(args_home+'config/pbmc.yaml')
args = namedtuple('Struct',params.keys())(*params.values())


if mode == 'nbr_net':

    # model_file = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+now.strftime('%Y%m%d%H%M')


    # logging.basicConfig(filename=model_file+'.log',
	# 					format='%(asctime)s %(levelname)-8s %(message)s',
	# 					level=logging.INFO,
	# 					datefmt='%Y-%m-%d %H:%M:%S')

    
    from plot import plt_umap 
    plt_umap.plot_umap_tcell(args)

    # from evals import results
    # results.topic_top_genes(args)

    # from evals import res_pbmcs
    # res_pbmcs.pbmc_sample_cells_with_latent(args)
    # res_pbmcs.topic_celltype_marker_genes(args)

elif mode == 'lr_net':

    # from model import interaction
    # interaction.run_model(args,model_file)

    # from model import interaction
    # interaction.load_model(args)

    # from model import interaction
    # interaction.eval_model(args)


    # from evals import results
    # results.topic_top_lr_genes(args)
    # results.topic_top_lr_pair_genes(args)

    from evals import res_pbmcs
    res_pbmcs.combine_topics_pbmc(args)

    # from evals import res_pbmcs
    # res_pbmcs.topic_celltype_marker_genes_lrnet(args)

