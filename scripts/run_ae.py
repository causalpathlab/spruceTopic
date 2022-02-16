import pandas as pd
import numpy as np
import preprocess 
import component_check 
from dlearn import ae
import datatable as dt
import importlib
import sys


def run_tcell_sc_data():

    data = "cd4_cd8_500cells_per_tissue_counts.txt.gz"
    # dtab = dt.fread("../output/tcell_data/"+data,sep="\t")
    # df = dtab.to_pandas()   
    df = pd.read_csv("../output/tcell_data/"+data,sep="\t")  

    df = df[ ["cell"]+\
            [x for x in df.columns if x not in["cell","sample"]]+\
            ["sample"]]
    df = preprocess.filter_minimal(df,cutoff=50)

    ##denoise tissue wise
    # df["sample"] = [x.split('_')[0] for x in df["sample"]]
    # tissues = [x.split('_')[0] for x in df["sample"].unique()]
    # df_xhat = pd.DataFrame()
    # for tissue in tissues:
    #     df_tissue = df[df["sample"]==tissue]
    #     df_tissue = df_tissue.reset_index(drop=True)
    #     y=df_tissue.iloc[:,-1]
    #     device = ae.torch.device('cuda' if ae.torch.cuda.is_available() else 'cpu')
    #     input_dims = df_tissue.iloc[:,1:-1].shape[1]
    #     print("Input dimension is "+ str(input_dims))

    #     batch_size = 64
    #     l_rate = 0.01
    #     epochs = 200
    #     layers = [128,32,16]
    #     latent_dims = 16

    #     model = ae.Autoencoder(input_dims,input_dims,latent_dims,layers).to(device)
    #     data = ae.load_data(df_tissue.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    #     print(model)
    #     loss_values = ae.train(model,data,device,epochs,l_rate)
    #     ae.get_encoded_h(df_tissue,model,device,"sc_epochs_"+str(epochs),loss_values)
    #     df_tissue_xhat = ae.get_xhat(df_tissue,model,device)
    #     df_xhat = df_xhat.append(df_tissue_xhat)
    #     del device
    #     del model
    #     del data
    #     print("completed..."+tissue)
    # df_xhat.to_csv("../output/tcell_data/cd4_cd8_500cells_per_tissue_counts_denoised.txt.gz",sep="\t",compression="gzip",index=False) 
    
    ##denoise all data
    y=df.iloc[:,-1]
    device = ae.torch.device('cuda' if ae.torch.cuda.is_available() else 'cpu')
    input_dims = df.iloc[:,1:-1].shape[1]
    print("Input dimension is "+ str(input_dims))

    batch_size = 64
    l_rate = 0.01
    epochs = 200
    layers = [1000,100]
    latent_dims = 100

    model = ae.Autoencoder(input_dims,input_dims,latent_dims,layers).to(device)
    data = ae.load_data(df.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    print(model)
    loss_values = ae.train(model,data,device,epochs,l_rate)
    ae.get_encoded_h(df,model,device,"sc_epochs_"+str(epochs),loss_values)
    df_xhat = ae.get_xhat(df,model,device)
    df_xhat.to_csv("../output/tcell_data/cd4_cd8_500cells_per_tissue_counts_denoised_allonce.txt.gz",sep="\t",compression="gzip",index=False) 

mode = sys.argv[1]
if mode =="tc_sc":
    run_tcell_sc_data()
