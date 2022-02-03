import pandas as pd
import numpy as np
from igraph import Graph
from igraph import plot as igplt


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



