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

    model_file = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+now.strftime('%Y%m%d%H%M')


    logging.basicConfig(filename=model_file+'.log',
					format='%(asctime)s %(levelname)-8s %(message)s',
					level=logging.INFO,
					datefmt='%Y-%m-%d %H:%M:%S')

    
    from plot import plt_umap 
    plt_umap.plot_umap_tcell(args)

    # from evals import results
    # results.topic_top_genes(args)

    # from evals import res_pbmcs
    # res_pbmcs.pbmc_sample_cells_with_latent(args)
    # res_pbmcs.topic_celltype_marker_genes(args)

elif mode == 'lr_net':

    # model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+now.strftime('%Y%m%d%H%M')
    # print('log - '+model_file)

    # logging.basicConfig(filename=model_file+'.log',
	# 				format='%(asctime)s %(levelname)-8s %(message)s',
	# 				level=logging.INFO,
	# 				datefmt='%Y-%m-%d %H:%M:%S')


    # from model import interaction
    # interaction.run_model(args,model_file)


    # from model import interaction
    # interaction.load_model(args)

    # from model import interaction
    # interaction.eval_model(args)

    # from evals import results
    # results.topic_top_lr_genes(args)
    # results.topic_top_lr_pair_genes(args)
    # results.assign_gene_bias(args)


    # from evals import res_tcells
    # res_tcells.combine_topics_tcells(args)


    from evals import res_pbmcs as ctr
    # from evals import res_tcells as ctr
    import igraph 

    g = ctr.cell_interaction_example(args)
    visual_style={}
    visual_style["vertex_label"] = g.vs["name"]
    visual_style["vertex_size"] = g.vs['size']
    visual_style["vertex_color"] = g.vs['color']
    visual_style["vertex_label_size"] = 8
    # visual_style["vertex_label_color"] = "darkblue"
    visual_style["edge_width"] = g.es['weight']
    visual_style["layout"] = g.layout_sugiyama()
    igraph.plot(g, target="network.pdf",**visual_style)
