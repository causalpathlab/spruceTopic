import pandas as pd
from gen_util.io import read_config


def filter_genes(df):
    
    # filter out the following-
    
    # 1) TODO non protein coding

    # 2)mitochondrial genes were eliminated previousl ?

    # 3) a gene was eliminated if the number of cells expressing this gene is <10.
    print("before filtering genes ",df.shape)
    drop_columns = [ col for col,val  in df.iloc[:,1:-1].sum(axis=0).iteritems() if val < 10 ]
    df = df.drop(drop_columns,axis=1)
    print("after filtering genes ",df.shape)

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
    
