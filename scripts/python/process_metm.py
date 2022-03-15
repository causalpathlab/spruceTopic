import pandas as pd
import numpy as np
from collections import namedtuple
from gen_util.io import read_config

config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/scmetm.yaml"
params = read_config(config)
args = namedtuple('Struct',params.keys())(*params.values())

def generate_gene_vals(df,top_n,top_genes,label):
    for x in range(df.shape[0]):
        gtab = df.T.iloc[:,x].sort_values(0,ascending=False)[:top_n].reset_index()
        gtab.columns = ["gene","val"]
        genes = gtab["gene"].values
        vals = gtab["val"].values
        for i in range(top_n):
            top_genes.append(["k"+str(x),label,"g"+str(i+1),genes[i],vals[i]])
    return top_genes


def get_top_genes(args,top_n=3):

    df_genes=pd.read_pickle(args.home+args.input+args.raw_genes)
    
    tag=args.model["data_tag"]
    immuneindex_file = args.home+args.input+args.model["sparse_label"]+"immune_index_"+tag+".npz"
    immunedat = np.load(immuneindex_file,allow_pickle=True)
    immuneindx = np.array([x[0] for x in immunedat["idx"]]).astype(np.int32)
    otherindx = np.array([x for x in range(df_genes.shape[0]) if x not in immuneindx]).astype(np.int32)

    df_beta1 = pd.read_csv(args.home+args.output+args.model["out"]+args.model["mfile"]+"etm_beta1_data.csv",sep="\t")
    df_beta2 = pd.read_csv(args.home+args.output+args.model["out"]+args.model["mfile"]+"etm_beta2_data.csv",sep="\t")

    df_beta1.columns = df_genes.iloc[immuneindx]["gene"].values
    df_beta2.columns = df_genes.iloc[otherindx]["gene"].values

    top_genes = []
    top_genes = generate_gene_vals(df_beta1,top_n,top_genes,"immune")
    top_genes = generate_gene_vals(df_beta2,top_n,top_genes,"non-immune")

    df_top_genes = pd.DataFrame(top_genes,columns=["Topic","GeneType","Genes","Gene","Proportion"])
    # top_gene_cols = ["t"+str(i) for i in range(top_n)]
    # top_gene_val_cols = ["v"+str(i) for i in range(top_n)]
    # df_top_genes[top_gene_cols] = df_top_genes["genes"].to_list()
    # df_top_genes[top_gene_val_cols] = df_top_genes["vals"].to_list()
    # df_top_genes = df_top_genes.drop(columns=["genes","vals"])
    df_top_genes.to_csv(args.home+args.output+args.model["mfile"]+"top_"+str(top_n)+"_genes_topic.tsv",sep="\t",index=False)

def sample_cells_with_latent(args,cell_n=500):

    from sklearn.cluster import KMeans
    df_h = pd.read_csv(args.home+args.output+args.model["out"]+args.model["mfile"]+"etm_hh_data.csv",sep="\t")
    kmeans = KMeans(n_clusters=df_h.shape[1]-1, random_state=0).fit(df_h.iloc[:,1:].to_numpy())
    df_h["cluster"] = kmeans.labels_
    df_h_sample = df_h.groupby("cluster").sample(cell_n, random_state=1)

    # df_h_sample = pd.DataFrame
    # for i in range(16):
    #     df_h = df_h.sort_values('hh'+str(i))
    #     if i ==0:
    #         df_h_sample = df_h.iloc[:10,:]
    #     else:
    #         df_h_sample = pd.concat([df_h_sample,df_h.iloc[:10,:]],axis=0, ignore_index=True)
    # df_h_sample = df_h_sample.drop_duplicates()
    df_h_sample.columns = [ x.replace("hh","k") for x in df_h_sample.columns]
    df_h_sample.to_csv(args.home+args.output+args.model["mfile"]+"hh_cell_topic_sample.tsv",sep="\t",index=False)

# get_top_genes(args)
sample_cells_with_latent(args)