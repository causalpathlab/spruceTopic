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

def plot_umap(args,etm='netm'):
	
	args_home = os.environ['args_home'] 
	dfz = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_'+etm+'_z.tsv.gz',sep='\t',compression='gzip')

	umap_2d = umap.UMAP(n_components=2, init='random', random_state=0,min_dist=0.1)
	proj_2d = umap_2d.fit(dfz.iloc[:,1:])
	dfz['umap1'] = umap_2d.embedding_[:,0]
	dfz['umap2'] = umap_2d.embedding_[:,1]
	
	dfh = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_'+etm+'_h.tsv.gz',sep='\t',compression='gzip')
	dfh.columns = [ x.replace('h','') for x in dfh.columns]
	dfh['topic'] = dfh.iloc[:,1:].idxmax(axis=1)
	dfz['Topic'] = pd.merge(dfz[['cell']],dfh[['cell','topic']],on='cell',how='left')['topic'].values

	
	cp = sns.color_palette(cc.glasbey, n_colors=len(dfz['Topic'].unique()))

	if len(dfz['Topic'].unique())<26:
		sns.scatterplot(data=dfz, x='umap1', y='umap2', hue='Topic',s=1,palette=cp)
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	else:
		sns.scatterplot(data=dfz, x='umap1', y='umap2', hue='Topic',s=1,palette=cp,legend=False)


	plt.savefig(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_'+etm+'_umap_z_label_topic.png',dpi=300)
	plt.close()


def plot_umap_with_annotation(args,label_index=1,etm='netm'):
	
	args_home = os.environ['args_home'] 
	dfz = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_'+etm+'_z.tsv.gz',sep='\t',compression='gzip')

	umap_2d = umap.UMAP(n_components=2, init='random', random_state=0,min_dist=0.1)
	proj_2d = umap_2d.fit(dfz.iloc[:,1:])
	dfz['umap1'] = umap_2d.embedding_[:,0]
	dfz['umap2'] = umap_2d.embedding_[:,1]

	dfl = pd.read_csv(args_home+args.input+args.sample+'_metadata.csv.gz',compression='gzip')
	label = dfl.columns[label_index]
	dfl = dfl.rename(columns={'Unnamed: 0':'cell'})
	dfz['Topic'] = pd.merge(dfz[['cell']],dfl[['cell',label]],on='cell',how='left')[label].values

	cp = sns.color_palette(cc.glasbey, n_colors=len(dfz['Topic'].unique()))
	t = sns.scatterplot(data=dfz, x='umap1', y='umap2', hue='Topic',s=1,palette=cp)
	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	t.axes.set_title("UMAP plot of cell-mixture",fontsize=30)
	t.set_xlabel("UMAP1",fontsize=20)
	t.set_ylabel("UMAP2",fontsize=20)
	plt.savefig(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_'+etm+'_umap_z_label_'+label+'.png',dpi=300)
	plt.close()

def plot_umap_with_annotation_mix(args,label_index=1,etm='netm'):
	
	args_home = os.environ['args_home'] 
	dfz = pd.read_csv(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_'+etm+'_z.tsv.gz',sep='\t',compression='gzip')

	umap_2d = umap.UMAP(n_components=2, init='random', random_state=0,min_dist=0.1)
	proj_2d = umap_2d.fit(dfz.iloc[:,1:])
	dfz['umap1'] = umap_2d.embedding_[:,0]
	dfz['umap2'] = umap_2d.embedding_[:,1]

	dfz['label'] = [x.split('_')[len(x.split('_'))-1] for x in dfz['cell']]

	f='/home/sishirsubedi/projects/spruce_topic/input/GSEmix/GSE176078_metadata.csv.gz'
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
	plt.savefig(args_home+args.output+args.nbr_model['out']+args.nbr_model['mfile']+'_'+etm+'_umap_z_label_'+label+'.png',dpi=300)
	plt.close()


