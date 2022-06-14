import os
from turtle import color
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [12.50, 10.50]
plt.rcParams['figure.autolayout'] = True
import umap
import umap.plot 
import igraph
import colorcet as cc
import seaborn as sns

def plot_umap(sp):
	
	dfz = sp.cell_topic.z

	umap_2d = umap.UMAP(n_components=2, init='random', random_state=0,min_dist=0.1)
	proj_2d = umap_2d.fit(dfz.iloc[:,1:])
	dfz['umap1'] = umap_2d.embedding_[:,0]
	dfz['umap2'] = umap_2d.embedding_[:,1]
	
	dfh = sp.cell_topic.h
	dfh.columns = [ x.replace('h','') for x in dfh.columns]
	dfh['topic'] = dfh.iloc[:,1:].idxmax(axis=1)
	dfz['Topic'] = pd.merge(dfz[['cell']],dfh[['cell','topic']],on='cell',how='left')['topic'].values

	
	cp = sns.color_palette(cc.glasbey, n_colors=len(dfz['Topic'].unique()))

	if len(dfz['Topic'].unique())<26:
		sns.scatterplot(data=dfz, x='umap1', y='umap2', hue='Topic',s=1,palette=cp)
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	else:
		sns.scatterplot(data=dfz, x='umap1', y='umap2', hue='Topic',s=1,palette=cp)


	plt.savefig(sp.model_id+'_umap_z_label_h_max.png',dpi=300)
	plt.close()


def plot_umap_with_annotation_mix(sp,label_index=8):
	
	dfz = sp.cell_topic.z

	umap_2d = umap.UMAP(n_components=2, init='random', random_state=0,min_dist=0.1)
	proj_2d = umap_2d.fit(dfz.iloc[:,1:])
	dfz['umap1'] = umap_2d.embedding_[:,0]
	dfz['umap2'] = umap_2d.embedding_[:,1]

	dfz['label'] = [x.split('_')[len(x.split('_'))-1] for x in dfz['cell']]

	f='/home/sishirsubedi/projects/data/GSE176078mix/GSE176078_metadata.csv.gz'
	dfl = pd.read_csv(f,compression='gzip')
	dfl = dfl.rename(columns={'Unnamed: 0':'cell'})

	dflabel = pd.DataFrame()
	dflabel['l1'] =  [x for x in dfz[dfz['label']=='GSE176078']['cell']]
	dflabel['l2'] =  [x.replace('_GSE176078','') for x in dfz[dfz['label']=='GSE176078']['cell']]
	dflabel = pd.merge(dflabel,dfl,right_on='cell',left_on='l2',how='left')

	label = dfl.columns[label_index]

	dfz = pd.merge(dfz,dflabel[['l1',label]],right_on='l1',left_on='cell',how='left')
	dfz[label] = dfz[label].mask(dfz[label].isna(), dfz['label'])
	
	cp = sns.color_palette(cc.glasbey, n_colors=len(dfz[label].unique()))
	t = sns.scatterplot(data=dfz, x='umap1', y='umap2', hue=label,s=1,palette=cp)
	plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
	# t.axes.set_title("UMAP plot of cell-mixture",fontsize=30)
	t.set_xlabel("UMAP1",fontsize=20)
	t.set_ylabel("UMAP2",fontsize=20)
	plt.savefig(sp.model_id+'_umap_z_label_'+label+'.png',dpi=300)
	plt.close()


