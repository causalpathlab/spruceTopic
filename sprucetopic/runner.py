import os 
import sys
import datetime
from _utils.io import read_config
from collections import namedtuple
import logging
import pandas as pd

mode= sys.argv[1]
now = datetime.datetime.now()
# args_home = '/home/BCCRC.CA/ssubedi/projects/spruce_topic/'
args_home = '/home/sishirsubedi/projects/spruce_topic/'

# os.chdir(args_home)
os.environ['args_home'] = args_home

params = read_config(args_home+'config/bcmix.yaml')
args = namedtuple('Struct',params.keys())(*params.values())

# import torch
# from scipy import sparse
# import numpy as np

# nbr_size = args.lr_model['train']['nbr_size']
# batch_size = args.lr_model['train']['batch_size']
# l_rate = args.lr_model['train']['l_rate']
# epochs = args.lr_model['train']['epochs']
# layers1 = args.lr_model['train']['layers1']
# layers2 = args.lr_model['train']['layers2']
# latent_dims = args.lr_model['train']['latent_dims']
# input_dims1 = args.lr_model['train']['input_dims1']
# input_dims2 = args.lr_model['train']['input_dims2']


# args_home = os.environ['args_home']
# device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# sparse_data = args_home+args.input+args.nbr_model['sparse_data']
# sparse_label = args_home+args.input+args.nbr_model['sparse_label']
# neighbour_data = args_home+ args.output+ args.lr_model['out']+args.nbr_model['mfile']+'_nbr.pkl'


# from model import interaction_topic as it
# import importlib
# importlib.reload(it)
# dl = it.load_data(sparse_data,sparse_label,neighbour_data, batch_size, device)
# train_dataloader =  dl.train_dataloader()
# x,y,z = next(iter(train_dataloader))


if mode == 'prep':

    model_file = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+now.strftime('%Y%m%d%H%M')


    logging.basicConfig(filename=model_file+'.log',
					format='%(asctime)s %(levelname)-8s %(message)s',
					level=logging.INFO,
					datefmt='%Y-%m-%d %H:%M:%S')
    
    from premodel import processing as prep     

    # prep.tenx_preprocessing(args)
    # prep.lr_preprocessing(args)
    
    prep.scanpy_processing(args)


elif mode == 'nbr_net':

    # model_file = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+now.strftime('%Y%m%d%H%M')


    # logging.basicConfig(filename=model_file+'.log',
	# 				format='%(asctime)s %(levelname)-8s %(message)s',
	# 				level=logging.INFO,
	# 				datefmt='%Y-%m-%d %H:%M:%S')

    
    # from model import netm
    # netm.run_model(args,model_file)

    from sprucetopic.model import cell_topic
    cell_topic.eval_model(args)

    from plot import plt_umap 
    plt_umap.plot_umap(args)
    plt_umap.plot_umap_with_annotation_mix(args,label_index=8)
    plt_umap.plot_umap_with_annotation_mix(args,label_index=1)
    plt_umap.plot_umap_with_annotation_mix(args,label_index=5)
    plt_umap.plot_umap_with_annotation_mix(args,label_index=6)
    plt_umap.plot_umap_with_annotation_mix(args,label_index=7)

    # from evals import results
    # results.topic_top_genes(args)
    # results.sample_cells_with_latent(args)


elif mode == 'lr_net':

    model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+now.strftime('%Y%m%d%H%M')
    print('log - '+model_file)

    logging.basicConfig(filename=model_file+'.log',
					format='%(asctime)s %(levelname)-8s %(message)s',
					level=logging.INFO,
					datefmt='%Y-%m-%d %H:%M:%S')

    # from model import ietm
    # ietm.generate_neighbours(args)

    from model import interaction_topic
    interaction_topic.run_model(args,model_file)

    # from model import ietm
    # ietm.load_model(args)

    # from model import ietm
    # ietm.eval_model(args)

    # from evals import results
    # results.topic_top_lr_genes(args)
    # results.topic_top_lr_pair_genes(args)
    # results.assign_gene_bias(args)
    # results.topics_summary(args)


    # from plot import plt_umap
    # plt_umap.plot_umap_pbmc(args)

