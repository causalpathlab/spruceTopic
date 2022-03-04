import pandas as pd
import numpy as np
import component_check 
from dlearn import metm
import importlib
import sys
import torch
from gen_util.io import read_config
import logging

def run_model():
    params = read_config()

    tag="10000"
    data_file = params["home"]+params["data"]+params["sparse_data"]+"_"+tag+".npz"
    meta_file = params["home"]+params["data"]+params["sparse_label"]+"_"+tag+".npz"
    immuneindex_file = params["home"]+params["data"]+params["sparse_label"]+"_immune_index_"+tag+".npz"

    logging.info("Starting...")
    logging.info(data_file)
    logging.info(meta_file)
    logging.info(immuneindex_file)
    #### etm tests

    batch_size = 64
    l_rate = 0.001
    epochs = 1000
    layers1 = [128,128,1000,32,16]
    layers2 = [128,128,1000,32,16]
    latent_dims = 16

    device = metm.torch.device('cuda' if metm.torch.cuda.is_available() else 'cpu')

    data = metm.load_sparse_data(data_file,meta_file,immuneindex_file, device,batch_size)

    input_dims1 = len(data.dataset.immuneindx)
    input_dims2 = data.dataset.shape[1] - input_dims1

    logging.info('Variables - \n batchsize:%s \n lrate:%s \n epochs:%s \n data:%s \n' % (str(batch_size), str(l_rate),str(epochs),str(data.dataset.shape)))
    logging.info("Input dimension - immune is "+ str(input_dims1))
    logging.info("Input dimension - others is "+ str(input_dims2))

    model = metm.ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
    logging.info(model)
    loss_values = metm.train(model,data,device,epochs,l_rate)

    fname = "_".join(str(x) for x in layers1) + "_" + str(epochs)+"_tcells_"+str(tag)
    torch.save(model.state_dict(), "../"+fname+"_etm.torch")
    dflv = pd.DataFrame(loss_values)
    dflv.to_csv("../"+fname+"_loss.txt",index=False)
    

def eval_model():
    params = read_config()

    tag="10000"
    data_file = params["home"]+params["data"]+params["sparse_data"]+"_"+tag+".npz"
    meta_file = params["home"]+params["data"]+params["sparse_label"]+"_"+tag+".npz"
    immuneindex_file = params["home"]+params["data"]+params["sparse_label"]+"_immune_index_"+tag+".npz"

    batch_size = 25000
    l_rate = 0.001
    epochs = 1000
    layers1 = [128,128,1000,32,16]
    layers2 = [128,128,1000,32,16]
    latent_dims = 16

    device = "cpu"
    fname = "_".join(str(x) for x in layers1) + "_" + str(epochs)+"_tcells_"+str(tag)

    data = metm.load_sparse_data(data_file,meta_file,immuneindex_file, device,batch_size)

    input_dims1 = len(data.dataset.immuneindx)
    input_dims2 = data.dataset.shape[1] - input_dims1

    model = metm.ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
    model.load_state_dict(torch.load("../"+fname+"_etm.torch"))
    model.eval()    

    dfloss = pd.read_csv("../"+fname+"_loss.txt")
    
    metm.get_encoded_h(data,model,fname,dfloss.iloc[:,0].values)

def check_z():
    runs = ["sc_raw_128_128_1000_32_16_1000_tcells_10000_data.csv"]
    for data in runs:
        df = pd.read_csv("../output/"+data,sep="\t")
        print(df.shape)
        # df = df.sample(n=25000)
        marker = ["CD4.c20.Treg.TNFRSF9",\
            "CD4.c06.Tm.ANXA1",\
            "CD4.c01.Tn.TCF7",\
            "CD8.c02.Tm.IL7R",\
            "CD8.c05.Tem.CXCR5",\
            "CD8.c07.Temra.CX3CR1"]

        y = np.array([x if x in marker else "others" for x in df.iloc[:,-1]])
        print(pd.Series(y).unique())
        component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"tsne",data.split(".")[0])
        component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"pca",data.split(".")[0])
mode = sys.argv[1]

if mode =="metm":
    logging.basicConfig(filename="../etm.log",
						format='%(asctime)s %(levelname)-8s %(message)s',
						level=logging.INFO,
						datefmt='%Y-%m-%d %H:%M:%S')
    run_model()
elif mode =="checkz":
    check_z()
elif mode =="metmeval":
    eval_model()

