import numpy as np
import pandas as pd
import annoy
import logging
logger = logging.getLogger(__name__)

class ApproxNN():
	def __init__(self, data, labels):
		self.dimension = data.shape[1]
		self.data = data.astype('float32')
		self.labels = labels

	def build(self, number_of_trees=50):
		self.index = annoy.AnnoyIndex(self.dimension,'angular')
		for i, vec in enumerate(self.data):
			self.index.add_item(i, vec.tolist())
		self.index.build(number_of_trees)

	def query(self, vector, k):
		indexes = self.index.get_nns_by_vector(vector.tolist(),k)
		return [self.labels[i][0] for i in indexes]

def get_NNmodels(df):
	model_list = {}
	for topic in df['topic'].unique():		
		model_ann = ApproxNN(df.loc[df['topic']==topic,list(df.columns[2:])].to_numpy(),df.loc[df['topic']==topic,['cell']].values )
		model_ann.build()
		model_list[topic] = model_ann
	return model_list

def get_neighbours(df,model_list,nbrsize,selected_topics):
		
	nbr_dict={}

	for idx,cell in enumerate(df['cell'].values):

		neighbours = np.array([ model_list[topic].query(df.iloc[idx,2:].values,k=nbrsize) for topic in selected_topics ]).flatten()
	
		## remove self 
		neighbours = neighbours[neighbours!=cell]

		if len(neighbours) != 159: 
			logger.info(str(idx)+'...'+str(len(neighbours)))
			neighbours = neighbours[1:]

		if len(neighbours) != 159: 
			logger.info('error...')
			logger.info(str(idx)+'...'+str(len(neighbours)))

		nbr_dict[cell] = df[df['cell'].isin(neighbours)].index.values

		if idx % 1000 == 0:
			logger.info(idx)

	return nbr_dict

def generate_neighbours(spruce,min_cells_in_topic=100):
	
	df_h = spruce.cell_topic.h

	df_l = spruce.data.raw_l_data
	df_r = spruce.data.raw_r_data
	df_l = df_l[df_l['index'].isin(df_h['cell'].values)]
	df_r = df_r[df_r['index'].isin(df_h['cell'].values)]

	dflatent = pd.merge(df_l['index'],df_h,how='left',left_on='index',right_on='cell')
	dflatent['cell'] = dflatent.iloc[:,2:].idxmax(axis=1)
	dflatent = dflatent.rename(columns={'cell':'topic','index':'cell'})

	selected_topics = dflatent.groupby('topic').filter( lambda x : len(x) >= min_cells_in_topic)['topic'].unique()

	model_list = get_NNmodels(dflatent.loc[dflatent['topic'].isin(selected_topics)])
	nbr_size = 5


	nbr_dict = get_neighbours(dflatent,model_list,nbr_size,selected_topics)
	return pd.DataFrame(nbr_dict).T

