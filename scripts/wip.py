import pandas as pd
import numpy as np
import preprocess 
import component_check 
from dlearn import etm
from dlearn import ae
import importlib


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
    epochs = 300
    layers = [9,6,3]
    latent_dims = 3

    model = etm.ETM(input_dims,input_dims,latent_dims,layers).to(device)
    data = etm.load_data(df.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    print(model)
    loss_values = etm.train(model,data,device,epochs,l_rate)
    etm.get_encoded_h(df,model,device,"wine_epochs_"+str(epochs),loss_values)

def run_toy_data():

    data = "toydata.csv"
    df = pd.read_csv("../output/toy_data/"+data,sep="\t")

    print("data shape before model", df.shape)
    y = np.array([1 if "CD4" in x else 0 for x in df.iloc[:,-1]])

    component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"pca")

    device = etm.torch.device('cuda' if etm.torch.cuda.is_available() else 'cpu')
    input_dims = df.iloc[:,1:-1].shape[1]
    print("Input dimension is "+ str(input_dims))

    batch_size = 32
    l_rate = 0.01
    epochs = 100
    layers = [10,2]
    latent_dims = 2

    model = etm.ETM(input_dims,input_dims,latent_dims,layers).to(device)
    data = etm.load_data(df.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    print(model)
    loss_values = etm.train(model,data,device,epochs,l_rate)
    etm.get_encoded_h(df,model,device,"toy_epochs_"+str(epochs),loss_values)



def run_sc_data():

    data = "cd4_cd8_500cells_per_tissue_counts.txt.gz"
    df = pd.read_csv("../output/cell_data/"+data,sep="\t")
    df = df[ ["cell"]+\
            [x for x in df.columns if x not in["cell","sample"]]+\
            ["sample"]]
    df = preprocess.filter_minimal(df)


    print("data shape before model", df.shape)
    y = np.array([1 if "CD4" in x else 0 for x in df.iloc[:,-1]])


    #### etm tests
    device = etm.torch.device('cuda' if etm.torch.cuda.is_available() else 'cpu')
    input_dims = df.iloc[:,1:-1].shape[1]
    print("Input dimension is "+ str(input_dims))

    batch_size = 64
    l_rate = 0.01
    epochs = 300
    layers = [128,32,16]
    latent_dims = 16

    model = etm.ETM(input_dims,input_dims,latent_dims,layers).to(device)
    data = etm.load_data(df.iloc[:,1:-1].to_numpy(),y,device,batch_size)
    print(model)
    loss_values = etm.train(model,data,device,epochs,l_rate)
    etm.get_encoded_h(df,model,device,"sc_epochs_"+str(epochs),loss_values)


run_sc_data()