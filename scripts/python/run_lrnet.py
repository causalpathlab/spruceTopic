from gen_util.io import read_config
from collections import namedtuple
from dlearn import lrnet
import pandas as pd
import logging

config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/scmetm.yaml"
params = read_config(config)
args = namedtuple('Struct',params.keys())(*params.values())
logging.basicConfig(filename=args.home+args.output+"lrnet.log",
						format='%(asctime)s %(levelname)-8s %(message)s',
						level=logging.INFO,
						datefmt='%Y-%m-%d %H:%M:%S')


nbr_size=101
batch_size=32
device = lrnet.torch.device('cuda' if lrnet.torch.cuda.is_available() else 'cpu')

data = lrnet.load_data(args,nbr_size,batch_size,device)


