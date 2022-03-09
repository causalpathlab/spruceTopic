from gen_util.io import read_config
from collections import namedtuple
import logging
import sys
from . import m_etm
config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/scmetm.yaml"
params = read_config(config)
args = namedtuple('Struct',params.keys())(*params.values())

logging.basicConfig(filename=args.model["mfile"]+".log",
						format='%(asctime)s %(levelname)-8s %(message)s',
						level=logging.INFO,
						datefmt='%Y-%m-%d %H:%M:%S')
 
mode = "train"
if mode == "train":
	m_etm.run_model(args)
elif mode == "eval":
	m_etm.eval_model(args)