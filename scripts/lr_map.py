import pandas as pd
import numpy as np
from igraph import Graph
from igraph import plot as igplt
from itertools import combinations



def plot_l_r_db():
    lr_db="../database/celltalkdb_v20220131_human_lr_pair.txt"
    df = pd.read_csv(lr_db,sep="\t")
    '''
    3398 entries in the db
    815 ligands
    780 receptors
    '''
    ### get gene names from the db
    genes = list(df.ligand_gene_symbol.unique())
    gene_type = ["blue" for x in genes]
    for receptor in df.receptor_gene_symbol.unique():
        if receptor not in genes:
            genes.append(receptor);gene_type.append("red")

    ## construct graph
    g = Graph()
    g.add_vertices(genes)
    g.vs['color'] = gene_type
    for i,row in df.iterrows():
        g.add_edge(row.ligand_gene_symbol,row.receptor_gene_symbol, weight=1)

    # del_nodes = [v.index for v in g.vs if v.degree()<5]
    visual_style={}
    visual_style["layout"] = g.layout_reingold_tilford()

    igplt(g, target='../output/interaction.pdf',**visual_style)



def plot_df(df):
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


