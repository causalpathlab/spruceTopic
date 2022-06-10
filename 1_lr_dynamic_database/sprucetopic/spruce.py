from dataloader import _sparse,_lr_augment,_lr_ddb
from model import _cell_topic,_interaction_topic,_interaction_topic_dir,_interaction_topic_Rdb
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import logging
logger = logging.getLogger(__name__)

class Spruce:
    def __init__(self):
        self.data = self.data()
        self.cell_topic = self.cell_topic()
        self.interaction_topic = self.interaction_topic()
        self.model_id = None
    
    class data:
        def __init__(self):
            self.raw_data = None 
            self.raw_data_genes = None 
            self.sparse_data = None 
            self.sparse_data_ids = None 

            self.raw_l_data = None
            self.raw_l_data_genes = None
            self.raw_r_data = None
            self.raw_r_data_genes = None
            self.raw_lr_data = None

            self.neighbour = None

    class cell_topic:
        def __init__(self):
            self.model = None
            self.z = None
            self.h = None
            self.beta = None

    class interaction_topic:
        def __init__(self):
            self.model = None
            self.z = None
            self.h = None
            self.beta1 = None
            self.beta2 = None

    def run_cell_topic(self,batch_size,l_rate,epochs,layers,latent_dims,device):

        logger.info('Starting model training...')

        data = _sparse.load_data(self.data.sparse_data,self.data.sparse_data_ids, device,batch_size)

        input_dims = data.dataset.shape[1]
        logging.info('Input dimension is '+ str(input_dims))

        self.cell_topic.model = _cell_topic.ETM(input_dims,latent_dims,layers).to(device)
        logging.info(self.cell_topic.model)
        loss_values = _cell_topic.train(self.cell_topic.model,data,epochs,l_rate)

        return loss_values

    def eval_cell_topic(self,batch_size,layers,latent_dims,device):

        logging.info('Starting model inference...')

        data = _sparse.load_data(self.data.sparse_data,self.data.sparse_data_ids, device,batch_size)

        input_dims = data.dataset.shape[1]
        logging.info('Input dimension is '+ str(input_dims))

        model = _cell_topic.ETM(input_dims,latent_dims,layers).to(device)
        model.load_state_dict(self.cell_topic.model)
        model.eval()

        return self.get_cell_topic_latent(data,model)

    def get_cell_topic_latent(self,data,model):

        for xx,y in data: break
        zz,m,v = model.encoder(xx)
        pr,hh = model.decoder(zz)
        hh = _cell_topic.torch.exp(hh)

        df_z = pd.DataFrame(zz.to('cpu').detach().numpy())
        df_z.columns = ['z'+str(i)for i in df_z.columns]
        df_z['cell'] = y
        df_z = df_z[ ['cell']+[x for x in df_z.columns if x not in['cell']]]

        df_h = pd.DataFrame(hh.to('cpu').detach().numpy())
        df_h.columns = ['h'+str(i)for i in df_h.columns]
        df_h['cell'] = y
        df_h = df_h[ ['cell']+[x for x in df_h.columns if x not in['cell']]]


        beta =  None
        for n,p in model.named_parameters():
            if n == 'decoder.lbeta':
                beta=p
        beta_smax = _cell_topic.nn.LogSoftmax(dim=-1)
        beta = _cell_topic.torch.exp(beta_smax(beta))

        df_beta = pd.DataFrame(beta.to('cpu').detach().numpy())

        return df_z,df_h,df_beta
    
    def run_interaction_topic(self,batch_size,epochs,layers1,layers2,latent_dims,input_dims1,input_dims2,device,f_loss):

        dl = _lr_ddb.load_data(self.data.cell_topic_h, self.data.raw_l_data, self.data.raw_r_data, self.data.neighbour, batch_size,device)

        train_dataloader =  dl.train_dataloader()

        logging.info('Input dimension - ligand is '+ str(input_dims1))
        logging.info('Input dimension - receptor is '+ str(input_dims2))
        model = _interaction_topic_Rdb.LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,f_loss)
        logging.info(model)

        trainer = _interaction_topic_Rdb.pl.Trainer(
        max_epochs=epochs,
        accelerator='gpu',
        plugins= _interaction_topic_Rdb.DDPPlugin(find_unused_parameters=False),
        gradient_clip_val=0.5,
        progress_bar_refresh_rate=50,
        enable_checkpointing=False)

        trainer.fit(model,train_dataloader)
        
        return model 

    def eval_interaction_topic(self,input_dims1,input_dims2,latent_dims,layers1,layers2):

        model = _interaction_topic_Rdb.LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,'temp.txt')

        model.load_state_dict(self.interaction_topic.model)
        model.eval()

        alpha,alpha_bias,beta,beta_bias =  None,None,None,None
        
        for n,p in model.named_parameters():
            print(n)
            if n == 'etm.decoder.p_alpha':
                alpha=p
            elif n == 'etm.decoder.p_beta':
                beta=p
            elif n == 'etm.decoder.alpha_bias':
                alpha_bias=p
            elif n == 'etm.decoder.beta_bias':
                beta_bias=p
		
        beta_smax = _interaction_topic_Rdb.nn.LogSoftmax(dim=-1)
        alpha = _interaction_topic_Rdb.torch.exp(beta_smax(alpha))
        beta = _interaction_topic_Rdb.torch.exp(beta_smax(beta))

        # df_alpha = pd.DataFrame(alpha.to('cpu').detach().numpy())
        # beta = beta.to('cpu').detach().numpy()
        # df_alpha_bias = pd.DataFrame(alpha_bias.to('cpu').detach().numpy())
        # beta_bias = beta_bias.to('cpu').detach().numpy()
        # return df_alpha,beta,df_alpha_bias,beta_bias

        df_beta = pd.DataFrame(beta.to('cpu').detach().numpy())
        alpha = alpha.to('cpu').detach().numpy()
        df_beta_bias = pd.DataFrame(beta_bias.to('cpu').detach().numpy())
        alpha_bias = alpha_bias.to('cpu').detach().numpy()
        return df_beta,alpha,df_beta_bias,alpha_bias


        
    
   
