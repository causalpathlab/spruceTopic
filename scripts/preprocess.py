import pandas as pd
import matplotlib.pylab as plt
from gen_util.io import read_config


def select_protein_coding_genes(df):
    params = read_config()
    pcg = params["home"]+params["database"]+params["protein_coding_genes"]
    df_pcg = pd.read_csv(pcg,header=None,names=["gene"])
    drop_columns = [x for x in df.columns[1:-1] if x not in df_pcg["gene"].values]
    df = df.drop(drop_columns,axis=1)
    return df

def scanpy_filter(df):
    
    import scanpy as sc
    
    ## convert df to anndata object
    obs = pd.DataFrame(df.iloc[:,0])
    var = pd.DataFrame(index=df.columns[1:-1])
    X = df.iloc[:,1:-1].to_numpy()
    adata = sc.AnnData(X, obs=obs, var=var, dtype='int32')

    ## filter scanpy default settings
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    ## scanpy pca check
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata, color='CST3')
    plt.savefig("../output/scanpy_test_pca.png");plt.close()

    ## get df back from anndata object
    df_scanpy = adata.to_df()
    df_scanpy.columns = adata.var_names
    df_scanpy["cell"] = adata.obs.cell
    dfjoin = pd.merge(df_scanpy,df[["cell","cluster"]],on="cell",how="left")
    dfjoin = dfjoin[["cell"]+[x for x in adata.var_names]+["cluster"]]
    return dfjoin


def filter_genes(df):
    
    print("initial data size : ",df.shape)

    # keeping protein coding genes
    df = select_protein_coding_genes(df)
    print("after selecting protein coding genes ",df.shape)

    # a gene was eliminated if the number of cells expressing this gene is <10.
    drop_columns = [ col for col,val  in df.iloc[:,1:-1].sum(axis=0).iteritems() if val < 10 ]
    df = df.drop(drop_columns,axis=1)
    print("after selecting genes with read counts >= 10 cells ",df.shape)
    
    return df


def read_sc_data():
    params = read_config()
    sample_path = params["home"]+params["data"]

    df_combine = pd.DataFrame()
    for sample in params["samples"]:
        print("processing--"+sample)
        df = pd.read_csv(sample_path+sample,sep="\t").T
        df = df.rename(columns=df.iloc[0])
        df = df.iloc[1:].reset_index()
        df = df.rename(columns={"index":"cell"})
        # df['sample'] = sample.split('_')[0]+"_"+sample.split('.')[1]
        if sample == params['sample_first']:
            df_combine = df
            print(df.shape,df_combine.shape)
        else:
            df_combine = pd.concat([df_combine, df], axis=0, ignore_index=True)
            print(df.shape,df_combine.shape)
    
    return df_combine


def select_cells():
    '''
    Create a test data for algorithm development from Zheng et al paper (GSE156728) 
    selected data - tissue type, cell type, cluster.
    '''
    
    celltype_cluster = ['CD4_c18','CD4_c13']
    tissue_type = ['BC', 'BCL', 'ESCA', 'MM', 'PACA', 'RC', 'THCA', 'UCEC']
    
    params = read_config()
    meta_path = params["home"]+params["database"]
    
    df = pd.read_csv(meta_path+params["metadata"],sep="\t")
    
    df["meta.cluster.sf"] =[x.split('.')[0]+"_"+x.split('.')[1] for x in df["meta.cluster"]]
    df = df[(df["meta.cluster.sf"].isin(celltype_cluster))]
    df = df[df["cancerType"].isin(tissue_type)]
    df['cluster'] = [x+"_"+y for x,y in zip(df["cancerType"],df["meta.cluster.sf"])]
    df = df[["cellID","cluster"]]
    df.to_csv("../output/all_tissue_cd4_c13_c18_cell_barcodes.csv",sep="\t",index=False)

def create_test_data():
    ## prepare data matrix for development
    df = read_sc_data()
    df = df.fillna(0)
    ## read selected cell barcodes
    df_bc = pd.read_csv("../output/all_tissue_cd4_c13_c18_cell_barcodes.csv",sep="\t")
    df_selected = pd.merge(df_bc,df,on="cell",how="left") 
    '''IMPORTANT:memory reaches 90%, 80k by 30k merge with 5k by 2'''
    df_selected= df_selected[ [ col for col in df_selected.columns if col != "cluster" ] + ["cluster"]]
    df_selected.to_csv("../output/cd4_c13_c18_counts.txt.gz",index=False,sep="\t",compression="gzip")
    
