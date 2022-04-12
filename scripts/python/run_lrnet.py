import sys
import os
from gen_util.io import read_config
from collections import namedtuple
from dlearn import lrnet
import pandas as pd
import logging
import datetime

now = datetime.datetime.now()
spath = os.path.dirname(__file__)
args_home = spath.replace('/scripts/python','/')

os.chdir(args_home)
os.environ['args_home'] = args_home

params = read_config(args_home+'/config/scmetm.yaml')
args = namedtuple('Struct',params.keys())(*params.values())
model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']+now.strftime('%Y%m%d%H%M')

print(model_file)
logging.basicConfig(filename=model_file+'.log',
						format='%(asctime)s %(levelname)-8s %(message)s',
						level=logging.INFO,
						datefmt='%Y-%m-%d %H:%M:%S')

def run_model(mode):
	nbr_size = args.lr_model['train']['nbr_size']
	batch_size = args.lr_model['train']['run_batch_size']
	l_rate = args.lr_model['train']['l_rate']
	epochs = args.lr_model['train']['epochs']
	layers1 = args.lr_model['train']['layers1']
	layers2 = args.lr_model['train']['layers2']
	latent_dims = args.lr_model['train']['latent_dims']

	if mode == 'generate_tensor':
		device = lrnet.torch.device('cuda' if lrnet.torch.cuda.is_available() else 'cpu')
		# lrnet.generate_tensors(args,nbr_size,device)
		lrnet.generate_tensors_nbrs_alltopic(args,nbr_size,device)

	elif mode=='train':
		device = lrnet.torch.device('cuda' if lrnet.torch.cuda.is_available() else 'cpu')
		data = lrnet.load_data(args,batch_size)

		input_dims1 = 539
		input_dims2 = 498
		logging.info('Input dimension - ligand is '+ str(input_dims1))
		logging.info('Input dimension - receptor is '+ str(input_dims2))

		model = lrnet.ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
		logging.info(model)
		loss_values = lrnet.train(model,data,epochs,l_rate)

		lrnet.torch.save(model.state_dict(), model_file+'etm.torch')
		dflv = pd.DataFrame(loss_values[0])
		dflv.to_csv(model_file+'loss.txt',index=False)
		dflv = pd.DataFrame(loss_values[1])
		dflv.to_csv(model_file+'loss2.txt',index=False)

	elif mode == 'eval':
		model_file = args_home+args.output+args.lr_model['out']+args.lr_model['mfile']
		device = 'cuda'
		batch_size = args.lr_model['eval']['batch_size']
		data = lrnet.load_data(args,batch_size)
		input_dims1 = 539
		input_dims2 = 498

		model = lrnet.ETM(input_dims1,input_dims2,latent_dims,layers1,layers2).to(device)
		model.load_state_dict(lrnet.torch.load(model_file+'etm.torch'))
		model.eval()	

		dfloss = pd.read_csv(model_file+'loss.txt',sep='\t')

		lrnet.get_latent(data,model,model_file,dfloss.iloc[:,0].values)

# def runner():
# 	mode='train'
# 	run_model(mode)

mode=sys.argv[1]
run_model(mode)
