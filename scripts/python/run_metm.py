import os 
import datetime
from gen_util.io import read_config
from collections import namedtuple
import logging
from dlearn import metm
import pandas as pd

now = datetime.datetime.now()
spath = os.path.dirname(__file__)
args_home = spath.replace('/scripts/python','/')

os.chdir(args_home)
os.environ['args_home'] = args_home

params = read_config(args_home+'/config/scmetm_v3.yaml')
args = namedtuple('Struct',params.keys())(*params.values())
model_file = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+now.strftime('%Y%m%d%H%M')

print(model_file)
logging.basicConfig(filename=model_file+'.log',
						format='%(asctime)s %(levelname)-8s %(message)s',
						level=logging.INFO,
						datefmt='%Y-%m-%d %H:%M:%S')

def run_model(args):

    logging.info('Starting model training...')
    
    data_file = args_home+args.input+args.nbr_model['sparse_data']
    meta_file = args_home+args.input+args.nbr_model['sparse_label']
    immuneindex_file = args_home+args.input+args.nbr_model['sparse_immune_label']
    
    batch_size = args.nbr_model['train']['batch_size']
    l_rate = args.nbr_model['train']['l_rate']
    epochs = args.nbr_model['train']['epochs']
    layers1 = args.nbr_model['train']['layers1']
    layers2 = args.nbr_model['train']['layers2']
    latent_dims = args.nbr_model['train']['latent_dims']
    
    # device = metm.torch.device('cuda' if metm.torch.cuda.is_available() else 'cpu')
    device = 'cuda:2'
    data = metm.load_sparse_data(data_file,meta_file,immuneindex_file, device,batch_size)

    input_dims1 = len(data.dataset.immuneindx)
    input_dims2 = data.dataset.shape[1] - input_dims1
    logging.info('Input dimension - immune is '+ str(input_dims1))
    logging.info('Input dimension - others is '+ str(input_dims2))

    model = metm.ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
    logging.info(model)
    loss_values = metm.train(model,data,epochs,l_rate)

    metm.torch.save(model.state_dict(), model_file+'etm.torch')
    dflv = pd.DataFrame(loss_values[0])
    dflv.to_csv(model_file+'loss.txt',index=False)
    dflv = pd.DataFrame(loss_values[1])
    dflv.to_csv(model_file+'loss2.txt',index=False)

def eval_model(args):
    
    logging.info('Starting model inference...')
    
    data_file = args_home+args.input+args.nbr_model['sparse_data']
    meta_file = args_home+args.input+args.nbr_model['sparse_label']
    immuneindex_file = args_home+args.input+args.nbr_model['sparse_immune_label']
    
    batch_size = args.nbr_model['eval']['batch_size']
    layers1 = args.nbr_model['train']['layers1']
    layers2 = args.nbr_model['train']['layers2']
    latent_dims = args.nbr_model['train']['latent_dims']
    
    model_file = args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']

    device = 'cpu'

    data = metm.load_sparse_data(data_file,meta_file,immuneindex_file, device,batch_size)

    print('data size is -',data.dataset.shape)

    input_dims1 = len(data.dataset.immuneindx)
    input_dims2 = data.dataset.shape[1] - input_dims1

    model = metm.ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
    model.load_state_dict(metm.torch.load(model_file+'etm.torch'))
    model.eval()    

    dfloss = pd.read_csv(model_file+'loss.txt')
    
    metm.get_latent(data,model,model_file,dfloss.iloc[:,0].values,'model')

run_model(args)