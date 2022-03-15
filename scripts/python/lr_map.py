import pandas as pd
import numpy as np
from igraph import Graph
from igraph import plot as igplt
import matplotlib.pylab as plt
from itertools import combinations

def genes_from_lr_db(mode,args):
    lr_db=args.home+args.database+args.lr_db
    df = pd.read_csv(lr_db,sep="\t")
    if mode == "gene":
        genes = []
        for gp in df["lr_pair"].values:
            gp = gp.split("_");genes.append(gp[0]);genes.append(gp[1])
        return genes
    elif mode == "lrpair":
        return df["lr_pair"].values

def keep_lr_genes(args):
    ### select only genes present in the ligand-receptor database 
    # need to run first time to filter out genes
    data = args.home+args.input+args.raw_data
    df_exp=pd.read_pickle(data)
    lr_genes = genes_from_lr_db("gene",args)
    selected_genes = [x for x in df_exp.columns if x in lr_genes ]
    df_exp = df_exp[["index"]+selected_genes]
    df_exp.to_pickle(args.home+args.input+args.raw_lr_data)

def graph_lr_db(data_genes,args):
    lr_db=args.home+args.database+args.lr_db
    df = pd.read_csv(lr_db,sep="\t")
    '''
    3398 entries in the db, 815 ligands, 780 receptors
    '''
    ### get gene names from the db
    genes = list(df.ligand_gene_symbol.unique())
    gene_type = ["DeepskyBlue" for x in genes]
    for receptor in df.receptor_gene_symbol.unique():
        if receptor not in genes:
            genes.append(receptor)
            gene_type.append("Medium Orchid")

    ## construct graph
    g = Graph()
    g.add_vertices(genes)
    g.vs['color'] = gene_type
    for i,row in df.iterrows():
        ligand = row.ligand_gene_symbol 
        receptor = row.receptor_gene_symbol 
        if ligand in data_genes and receptor in data_genes:
            g.add_edge(row.ligand_gene_symbol,row.receptor_gene_symbol, weight=0)

    del_nodes = [v.index for v in g.vs if v.degree()<1]
    g.delete_vertices(del_nodes)
    return g

def graph_df(df):
    genes =  [x for x in df.T.index[0:-1]]
    ## construct graph
    g = Graph()
    g.add_vertices(genes)
    gedges =[]
    for i,row in df.iloc[0:-1,:-1].iterrows():
        pair_indexes = [comb for comb in combinations(row[row==1].index, 2)]
        for pair_index in pair_indexes:
            if pair_index not in gedges:
                g.add_edge(pair_index[0],pair_index[1],weight=1)
                gedges.append(pair_index)
            elif pair_index in gedges:
                eid = g.get_eid(pair_index[0],pair_index[1])
                g.es[eid]["weight"] += 1
    return g

def update_graph(g,df_mat):
    all_edges = {}
    for e in g.es:
        all_edges[g.vs[e.source]['name']+'_'+g.vs[e.target]['name']]=e.index
    del_edges = [ all_edges[x] for x in all_edges.keys() if x not in df_mat.index ]
    g.delete_edges(del_edges)
    del_nodes = [v.index for v in g.vs if v.degree()<1]
    g.delete_vertices(del_nodes)

    ## update weight
    for i in df_mat.index:
        s = i.split('_')[0]
        t = i.split('_')[1]
        eid = g.get_eid(s,t)
        for c in df_mat.columns:
            # if df_mat.loc[i,c] > 1:
            #     g.es[eid]["weight"] += 1
            g.es[eid]["weight"] += df_mat.loc[i,c]/df_mat.loc[:,c].sum()

    return g

def plot_graph(g,df,pltname):
    clustn = len(df["cluster"].unique())
    colors = [plt.cm.Dark2(x) for x in range(clustn)]
    color_dict = dict(zip(range(clustn), colors))

    edge_colors = []
    for e in g.es:
        id =g.vs[e.source]['name']+'_'+g.vs[e.target]['name']
        edge_colors.append(df.loc[id,'cluster'])
        

    edge_colors = [color_dict[int(x)]for x in edge_colors]
    g.es['color'] = edge_colors
    visual_style={}
    visual_style["vertex_label"] = g.vs["name"]
    visual_style["vertex_size"] = 10
    visual_style["vertex_label_size"] = 8
    visual_style["vertex_label_color"] = "black"
    visual_style["edge_width"] = g.es["weight"]
    visual_style["layout"] = g.layout_circle()
    visual_style["vertex_label_dist"] = 1
    visual_style["vertex_label_angle"] = 45    
    igplt(g, target=pltname+"_network.pdf",**visual_style)
