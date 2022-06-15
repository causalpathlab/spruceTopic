from dataloader import _sparse,_lr_augment
from model import _cell_topic,_interaction_topic
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
            self.beta_mean = None
            self.beta_var = None

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
        loss_values = _cell_topic.train(self.cell_topic.model,data,epochs,l_rate,batch_size)

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
        bm,bv,theta,beta  = model.decoder(zz)

        df_z = pd.DataFrame(zz.to('cpu').detach().numpy())
        df_z.columns = ['z'+str(i)for i in df_z.columns]
        df_z['cell'] = y
        df_z = df_z[ ['cell']+[x for x in df_z.columns if x not in['cell']]]

        df_h = pd.DataFrame(theta.to('cpu').detach().numpy())
        df_h.columns = ['h'+str(i)for i in df_h.columns]
        df_h['cell'] = y
        df_h = df_h[ ['cell']+[x for x in df_h.columns if x not in['cell']]]


        beta_mean =  None
        beta_var = None
        for n,p in model.named_parameters():
            if n == 'decoder.beta_mean':
                beta_mean = p
            if n == 'decoder.beta_lnvar':
                beta_var = p

        df_beta = pd.DataFrame(beta_mean.to('cpu').detach().numpy())
        df_beta_var = pd.DataFrame(beta_var.to('cpu').detach().numpy())

        return df_z,df_h,df_beta,df_beta_var
    
    def run_interaction_topic(self,batch_size,epochs,layers1,layers2,latent_dims,input_dims1,input_dims2,device,f_loss):

        dl = _lr_augment.load_data(self.cell_topic.h, self.data.raw_l_data, self.data.raw_r_data, self.data.raw_lr_data, self.data.neighbour, batch_size,device)

        train_dataloader =  dl.train_dataloader()

        logging.info('Input dimension - ligand is '+ str(input_dims1))
        logging.info('Input dimension - receptor is '+ str(input_dims2))
        model = _interaction_topic.LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,f_loss)
        logging.info(model)

        trainer = _interaction_topic.pl.Trainer(
        max_epochs=epochs,
        accelerator='gpu',
        plugins= _interaction_topic.DDPPlugin(find_unused_parameters=False),
        gradient_clip_val=0.5,
        progress_bar_refresh_rate=50,
        enable_checkpointing=False)

        trainer.fit(model,train_dataloader)
        
        return model 

    def eval_interaction_topic(self,input_dims1,input_dims2,latent_dims,layers1,layers2):

        model = _interaction_topic.LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,'temp.txt')

        model.load_state_dict(self.interaction_topic.model)
        model.eval()

        beta1,beta1_bias,beta2,beta2_bias =  None,None,None,None
        
        for n,p in model.named_parameters():
            print(n)
            if n == 'etm.decoder.lbeta1':
                beta1=p
            elif n == 'etm.decoder.lbeta2':
                beta2=p
            elif n == 'etm.decoder.lbeta1_bias':
                beta1_bias=p
            elif n == 'etm.decoder.lbeta2_bias':
                beta2_bias=p
            

        beta_smax = _interaction_topic.nn.LogSoftmax(dim=-1)
        beta1 = _interaction_topic.torch.exp(beta_smax(beta1))
        beta2 = _interaction_topic.torch.exp(beta_smax(beta2))

        df_beta1 = pd.DataFrame(beta1.to('cpu').detach().numpy())
        df_beta2 = pd.DataFrame(beta2.to('cpu').detach().numpy())
        df_beta1_bias = pd.DataFrame(beta1_bias.to('cpu').detach().numpy())
        df_beta2_bias = pd.DataFrame(beta2_bias.to('cpu').detach().numpy())

        return df_beta1,df_beta2,df_beta1_bias,df_beta2_bias

    def interactions_prob(self,input_dims1,input_dims2,latent_dims,layers1,layers2):

        model = _interaction_topic.LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,'_txt_')
        model.load_state_dict(self.interaction_topic.model)
        model.eval()

        df_h = self.cell_topic.h

        df_l = self.data.raw_l_data
        df_r = self.data.raw_r_data
        df_l = df_l[df_l['index'].isin(df_h['cell'].values)]
        df_r = df_r[df_r['index'].isin(df_h['cell'].values)]
        df_lr = self.data.raw_lr_data
        df_lr = df_lr.loc[df_l.columns[1:],df_r.columns[1:]]
        df_nbr = self.data.neighbour

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        nbrmat = torch.tensor(df_nbr.values.astype(np.compat.long),requires_grad=False).to(device)
        lmat = torch.tensor(df_l.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(device)
        rmat = torch.tensor(df_r.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(device)
        lrmat = torch.tensor(df_lr.values.astype(np.float32),requires_grad=False).to(device)

        topics_prob = []
        for idx in range(df_nbr.shape[0]):

            cm_l = lmat[idx].unsqueeze(0)
            cm_r = rmat[idx].unsqueeze(0)

            cm_lr = torch.mm(cm_l,lrmat).mul(cm_r)
            cm_rl = torch.mm(cm_r,torch.t(lrmat)).mul(cm_l)

            nbr_idxs = nbrmat[idx]
            
            cn_l = lmat[nbr_idxs]
            cn_r = rmat[nbr_idxs]

            lprime,rprime =  cm_lr + torch.mm(cn_l,lrmat).mul(cn_r) , cm_rl + torch.mm(cn_r,torch.t(lrmat)).mul(cn_l)

            pr,m1,v1,m2,v2,h = model.etm(lprime,rprime)

            h_smax = nn.LogSoftmax(dim=-1)
            h = torch.exp(h_smax(h))

            topics_prob.append(h.detach().numpy())
        
        return topics_prob

    def interactions_state_summary(self,input_dims1,input_dims2,latent_dims,layers1,layers2):

        model = _interaction_topic.LitETM(input_dims1,input_dims2,latent_dims,layers1,layers2,'_txt_')
        model.load_state_dict(self.interaction_topic.model)
        model.eval()

        df_h = self.cell_topic.h

        df_l = self.data.raw_l_data
        df_r = self.data.raw_r_data
        df_l = df_l[df_l['index'].isin(df_h['cell'].values)]
        df_r = df_r[df_r['index'].isin(df_h['cell'].values)]
        df_lr = self.data.raw_lr_data
        df_lr = df_lr.loc[df_l.columns[1:],df_r.columns[1:]]
        df_nbr = self.data.neighbour

        # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        device='cpu'
        nbrmat = torch.tensor(df_nbr.values.astype(np.compat.long),requires_grad=False).to(device)
        lmat = torch.tensor(df_l.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(device)
        rmat = torch.tensor(df_r.iloc[:,1:].values.astype(np.float32),requires_grad=False).to(device)
        lrmat = torch.tensor(df_lr.values.astype(np.float32),requires_grad=False).to(device)

        topics = []
        for idx in range(df_nbr.shape[0]):

            cm_l = lmat[idx].unsqueeze(0)
            cm_r = rmat[idx].unsqueeze(0)

            cm_lr = torch.mm(cm_l,lrmat).mul(cm_r)
            cm_rl = torch.mm(cm_r,torch.t(lrmat)).mul(cm_l)

            nbr_idxs = nbrmat[idx]
            
            cn_l = lmat[nbr_idxs]
            cn_r = rmat[nbr_idxs]

            lprime,rprime =  cm_lr + torch.mm(cn_l,lrmat).mul(cn_r) , cm_rl + torch.mm(cn_r,torch.t(lrmat)).mul(cn_l)

            pr,m1,v1,m2,v2,h = model.etm(lprime,rprime)

            h_smax = nn.LogSoftmax(dim=-1)
            h = torch.exp(h_smax(h))

            topics.append(list(pd.DataFrame(h.detach().numpy()).idxmax(axis=1).values))
        
        df_it = pd.DataFrame(topics)
        df_it['cell'] = df_l['index']
        df_it = df_it[['cell']+[x for x in df_it.columns[:-1]]]
        df_it.to_csv(self.model_id+'_ietm_interaction_states.tsv.gz',sep='\t',index=False,compression='gzip')
            
