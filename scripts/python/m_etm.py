import pandas as pd
import numpy as np
from scmetm import component_check 
from dlearn import  metm
import torch
import logging
logger = logging.getLogger(__name__)

def run_model(args):

    logger.info("Starting model training...")
    
    tag=args.model["data_tag"]
    data_file = args.home+args.input+args.model["sparse_data"]+tag+".npz"
    meta_file = args.home+args.input+args.model["sparse_label"]+tag+".npz"
    immuneindex_file = args.home+args.input+args.model["sparse_label"]+"immune_index_"+tag+".npz"
    
    batch_size = args.model["train"]["batch_size"]
    l_rate = args.model["train"]["l_rate"]
    epochs = args.model["train"]["epochs"]
    layers1 = args.model["train"]["layers1"]
    layers2 = args.model["train"]["layers2"]
    latent_dims = args.model["train"]["latent_dims"]
    
    model_file = args.home+args.output+args.model["out"]+args.model["mfile"]

    device = metm.torch.device('cuda' if metm.torch.cuda.is_available() else 'cpu')

    data = metm.load_sparse_data(data_file,meta_file,immuneindex_file, device,batch_size)

    input_dims1 = len(data.dataset.immuneindx)
    input_dims2 = data.dataset.shape[1] - input_dims1
    logger.info("Input dimension - immune is "+ str(input_dims1))
    logger.info("Input dimension - others is "+ str(input_dims2))

    model = metm.ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
    logger.info(model)
    loss_values = metm.train(model,data,device,epochs,l_rate)

    torch.save(model.state_dict(), model_file+"etm.torch")
    dflv = pd.DataFrame(loss_values)
    dflv.to_csv(model_file+"loss.txt",index=False)
    

def eval_model(args):
    
    logger.info("Starting model inference...")
    
    tag=args.model["data_tag"]
    data_file = args.home+args.input+args.model["sparse_data"]+tag+".npz"
    meta_file = args.home+args.input+args.model["sparse_label"]+tag+".npz"
    immuneindex_file = args.home+args.input+args.model["sparse_label"]+"immune_index_"+tag+".npz"
    
    batch_size = args.model["eval"]["batch_size"]
    l_rate = args.model["eval"]["l_rate"]
    epochs = args.model["eval"]["epochs"]
    layers1 = args.model["eval"]["layers1"]
    layers2 = args.model["eval"]["layers2"]
    latent_dims = args.model["eval"]["latent_dims"]
    
    model_file = args.home+args.output+args.model["out"]+args.model["mfile"]

    device = "cpu"

    data = metm.load_sparse_data(data_file,meta_file,immuneindex_file, device,batch_size)

    print("data size is -",data.dataset.shape)

    input_dims1 = len(data.dataset.immuneindx)
    input_dims2 = data.dataset.shape[1] - input_dims1

    model = metm.ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
    model.load_state_dict(torch.load(model_file+"etm.torch"))
    model.eval()    

    dfloss = pd.read_csv(model_file+"loss.txt")
    
    metm.get_latent(data,model,model_file,dfloss.iloc[:,0].values,"model")

def check_z(args):
    model_file = args.home+args.output+args.model["out"]+args.model["mfile"]

    df = pd.read_csv(model_file+"data.csv",sep="\t")
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
    component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"tsne",model_file)
    component_check.run_component_analysis(df.iloc[:,1:df.shape[1]-1],y,"pca",model_file)


def check_z_umap(args):

    import scanpy as sc
    import numpy as np
    import matplotlib.pylab as plt
    plt.rcParams["figure.figsize"] = [12.50, 10.50]
    plt.rcParams["figure.autolayout"] = True

    model_file = args.home+args.output+args.model["out"]+args.model["mfile"]
    df = pd.read_csv(model_file+"data.csv",sep="\t")
    print(df.shape)
    marker = ["CD4.c20.Treg.TNFRSF9",\
        "CD4.c06.Tm.ANXA1",\
        "CD4.c01.Tn.TCF7",\
        "CD8.c02.Tm.IL7R",\
        "CD8.c05.Tem.CXCR5",\
        "CD8.c07.Temra.CX3CR1"]

    df["sample"] = np.array([x if x in marker else "others" for x in df.iloc[:,-1]])

    obs = pd.DataFrame(df.iloc[:,0])
    var = pd.DataFrame(index=df.columns[1:-1])
    X = df.iloc[:,1:-1].to_numpy()
    adata = sc.AnnData(X, obs=obs, var=var, dtype='int32')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    adata.obs["pclust_mod"]=df["sample"].values
    sc.pl.umap(adata, color=['pclust_mod'],s=10)
    plt.savefig(model_file+"etm_umap.png");plt.close()


