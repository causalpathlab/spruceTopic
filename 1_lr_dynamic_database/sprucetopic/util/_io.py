
def read_config(config_file):
	import yaml
	with open(config_file) as f:
		params = yaml.safe_load(f)
	return params

def pickle_obj(f_name, data):
	import pickle
	pikd = open(f_name + '.pickle', 'wb')
	pickle.dump(data, pikd)
	pikd.close()

def unpickle_obj(f_name):
	import pickle
	pikd = open(f_name, 'rb')
	data = pickle.load(pikd)
	pikd.close()
	return data

def compress_pickle_obj(f_name, data):
	import pickle 
	import bz2
	with bz2.BZ2File(f_name + '.pbz2', 'w') as f:
		pickle.dump(data, f)

def decompress_pickle_obj(f_name):
	import pickle 
	import bz2
	data = bz2.BZ2File(f_name, 'rb')
	data = pickle.load(data)
	return data