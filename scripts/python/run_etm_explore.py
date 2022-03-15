import pandas as pd
import numpy as np
import preprocess 
import component_check 
from dlearn import etm
import datatable as dt
import sys


def run_wine_data():

    data = "wine_data.csv"
    df = pd.read_csv("../output/wine_data/"+data,sep="\t")

    drop_columns_sample = [ col for col,val  in df.iloc[:,1:-1].sum(axis=0).iteritems() if val < 1 ]
    df = df.drop(drop_columns_sample,axis=1)

    print("data shape before model", df.shape)
    y = np.array([1 if "CD4" in x else 0 for x in df.iloc[:,-1]])
    device = etm.torch.device('cuda' if etm.torch.cuda.is_available() else 'cpu')
    input_dims = df.iloc[:,1:-1].shape[1]
    print("Input dimension is "+ str(input_dims))

    batch_size = 32
    l_rate = 0.01
    epochs = 1000
    layers = [9,6,3]
    latent_dims = 3

    model = etm.ETM(input_dims,input_dims,latent_dims,layers).to(device)
    data = etm.load_data(df.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    print(model)
    loss_values = etm.train(model,data,device,epochs,l_rate)
    etm.get_encoded_h(df,model,device,"wine_epochs_"+str(epochs),loss_values)

def run_toy_data():

    data = "toydatav2.csv"
    df = pd.read_csv("../output/toy_data/"+data,sep="\t")

    print("data shape before model", df.shape)
    y = np.array([1 if "CD4" in x else 0 for x in df.iloc[:,-1]])

    component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"pca")

    device = etm.torch.device('cuda' if etm.torch.cuda.is_available() else 'cpu')
    input_dims = df.iloc[:,1:-1].shape[1]
    print("Input dimension is "+ str(input_dims))

    batch_size = 32
    l_rate = 0.01
    epochs = 1000
    layers = [25,3]
    latent_dims = 3

    model = etm.ETM(input_dims,input_dims,latent_dims,layers).to(device)
    data = etm.load_data(df.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    print(model)
    loss_values = etm.train(model,data,device,epochs,l_rate)
    etm.get_encoded_h(df,model,device,"toy_epochs_"+str(epochs),loss_values)

def run_tcell_sc_data():

    data = "cd4_cd8_500cells_per_tissue_counts.txt.gz"
    df = pd.read_csv("../output/tcell_data/"+data,sep="\t")
    df = df[ ["cell"]+\
            [x for x in df.columns if x not in["cell","sample"]]+\
            ["sample"]]
    # df = preprocess.filter_minimal(df)

    print("data shape before model", df.shape)
    y = np.array([1 if "CD4" in x else 0 for x in df.iloc[:,-1]])

    y=df.iloc[:,-1]
    # component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"pca")

    #### etm tests
    device = etm.torch.device('cuda' if etm.torch.cuda.is_available() else 'cpu')
    input_dims = df.iloc[:,1:-1].shape[1]
    print("Input dimension is "+ str(input_dims))

    batch_size = 64
    l_rate = 0.001
    epochs = 2000
    layers = [128,128,32,16]
    latent_dims = 16

    model = etm.ETM(input_dims,input_dims,latent_dims,layers).to(device)
    data = etm.load_data(df.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    print(model)
    loss_values = etm.train(model,data,device,epochs,l_rate)
    df_z = etm.get_encoded_h(df,model,device,"tcell_epochs_"+str(epochs),loss_values)

    component_check.run_component_analysis(df_z.iloc[:,1:df_z.shape[1]-1],df_z.iloc[:,-1],"pca")
    component_check.run_component_analysis(df_z.iloc[:,1:df_z.shape[1]-1],df_z.iloc[:,-1],"tsne")

def run_pbmc_sc_data():

    data = "pbmc3k_scanpy_raw_counts.txt.gz"

    dtab = dt.fread("../output/pbmc_data/"+data)
    df = dtab.to_pandas()   
    
    df = df.rename(columns={"index":"cell"})
    df['sample'] = "pbmc"


    df = preprocess.filter_minimal(df,50)
    print("data shape before model", df.shape)

    dfmeta = pd.read_csv("../output/pbmc_data/pbmc_metadata.csv",sep="\t")
    df = pd.merge(dfmeta,df,on="cell",how="left")
    df["sample_y"] = df["sample_x"]
    df = df.drop("sample_x",axis=1)
    df = df.rename(columns={"sample_y":"sample"})


    y=df.iloc[:,-1]
    # component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"pca")

    #### etm tests
    device = etm.torch.device('cuda' if etm.torch.cuda.is_available() else 'cpu')
    input_dims = df.iloc[:,1:-1].shape[1]
    print("Input dimension is "+ str(input_dims))

    batch_size = 64
    l_rate = 0.01
    epochs = 2000
    layers = [128,64,10]
    latent_dims = 10

    model = etm.ETM(input_dims,input_dims,latent_dims,layers).to(device)
    data = etm.load_data(df.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    print(model)
    loss_values = etm.train(model,data,device,epochs,l_rate)

    df_z = etm.get_encoded_h(df,model,device,"pbmc_epochs_"+str(epochs),loss_values)

    component_check.run_component_analysis(df_z.iloc[:,1:df_z.shape[1]-1],df_z.iloc[:,-1],"pca")
    component_check.run_component_analysis(df_z.iloc[:,1:df_z.shape[1]-1],df_z.iloc[:,-1],"tsne")

def run_bc_sc_data():

    data = "GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"
    dtab = dt.fread("../input/bc_data/"+data)
    df = dtab.to_pandas()   

    select_cells = [ x for x in df.columns if "Pooled" not in x]
    df = df[df['gene_type']=="protein_coding"]
    df = df[select_cells]
    df = df.drop(["gene_id","gene_type"],axis=1)
    df = df.T
    df = df.rename(columns=df.iloc[0])
    df = df.iloc[1:].reset_index()
    df = df.rename(columns={"index":"cell"})
    df['sample'] = [x.split('_')[0] for x in df["cell"]]

    ## save this as filtered count data 
    # df.to_csv("../output/bc_sc_data/bc_counts.txt.gz",index=False,sep="\t",compression="gzip")

    df = preprocess.filter_minimal(df,cutoff=50)
    print("data shape before model", df.shape)

    y=df.iloc[:,-1]
    # component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"pca")


    #### etm tests
    device = etm.torch.device('cuda' if etm.torch.cuda.is_available() else 'cpu')
    input_dims = df.iloc[:,1:-1].shape[1]
    print("Input dimension is "+ str(input_dims))

    batch_size = 64
    l_rate = 0.01
    epochs = 2000
    layers = [128,50,10]
    latent_dims = 10

    model = etm.ETM(input_dims,input_dims,latent_dims,layers).to(device)
    data = etm.load_data(df.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    print(model)
    loss_values = etm.train(model,data,device,epochs,l_rate)

    df_z = etm.get_encoded_h(df,model,device,"bc_epochs_"+str(epochs),loss_values)

    component_check.run_component_analysis(df_z.iloc[:,1:df_z.shape[1]-1],df_z.iloc[:,-1],"pca","fpca")
    component_check.run_component_analysis(df_z.iloc[:,1:df_z.shape[1]-1],df_z.iloc[:,-1],"tsne","ftsne")

def run_nyt_data():

    data = "nydata.csv"
    df = pd.read_csv("../output/news_data/"+data,sep="\t")

    print("data shape before model", df.shape)
    y = np.array([1 if "CD4" in x else 0 for x in df.iloc[:,-1]])

    # component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"pca")

    device = etm.torch.device('cuda' if etm.torch.cuda.is_available() else 'cpu')
    input_dims = df.iloc[:,1:-1].shape[1]
    print("Input dimension is "+ str(input_dims))

    batch_size = 64
    l_rate = 0.01
    epochs = 2000
    layers = [128,32,16]
    latent_dims = 16

    model = etm.ETM(input_dims,input_dims,latent_dims,layers).to(device)
    data = etm.load_data(df.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    print(model)
    loss_values = etm.train(model,data,device,epochs,l_rate)
    etm.get_encoded_h(df,model,device,"nyt_epochs_"+str(epochs),loss_values)

mode = sys.argv[1]
if mode =="toy":
    run_toy_data()
elif mode =="tc_sc":
    run_tcell_sc_data()
elif mode =="pbmc_sc":
    run_pbmc_sc_data()
elif mode =="bc_sc":
    run_bc_sc_data()
elif mode == "wine":
    run_wine_data()
elif mode == "nyt":
    run_nyt_data()