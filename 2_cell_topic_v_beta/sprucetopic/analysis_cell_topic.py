

import experiment_cell_topic 
import pandas as pd


experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'
spruce = experiment_cell_topic.get_experiment(experiment_home)
print(spruce.model_id)

##################################################################


##################################################################
                 ## umap figures ##
##################################################################

# generate umap coordinates based on softmaxed latent dimensions 
import umap

dfh = spruce.cell_topic.h.copy()
df_umap= pd.DataFrame()
df_umap['cell'] = dfh['cell']

umap_2d = umap.UMAP(n_components=2, init='random', random_state=0,min_dist=0.0,metric='cosine')
proj_2d = umap_2d.fit(dfh.iloc[:,1:])

df_umap[['umap1','umap2']] = umap_2d.embedding_[:,[0,1]]
df_umap['topic'] = [x.replace('h','') for x in dfh.iloc[:,1:].idxmax(axis=1)]
df_umap.to_csv(spruce.model_id+'_ct_h_umap_cordinates.csv.gz',index=False,compression='gzip')

##################################################################

import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [15, 10]
plt.rcParams['figure.autolayout'] = True
import colorcet as cc
import seaborn as sns

# get previously saved umap coordinates

df_umap = pd.read_csv(spruce.model_id+'_ct_h_umap_cordinates.csv.gz',compression='gzip')

##################################################################

# plot umap with argmax topic

cp = sns.color_palette(cc.glasbey, n_colors=len(df_umap['topic'].unique()))
p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='topic',s=1,palette=cp,legend=False)
# plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
# t.axes.set_title("UMAP plot of cell-mixture",fontsize=30)
p.set_xlabel("UMAP1",fontsize=20)
p.set_ylabel("UMAP2",fontsize=20)
plt.savefig(spruce.model_id+'_ct_umap_h_argmax_label.png');plt.close()

##################################################################

#plot umap with cell type major label

df_umap['label'] = [x.split('_')[len(x.split('_'))-1] for x in df_umap['cell']]

f='/home/sishirsubedi/projects/data/GSE176078mix/GSE176078_metadata.csv.gz'
dfl = pd.read_csv(f,compression='gzip')
dfl = dfl.rename(columns={'Unnamed: 0':'cell'})

dflabel = pd.DataFrame()
dflabel['l1'] =  [x for x in df_umap[df_umap['label']=='GSE176078']['cell']]
dflabel['l2'] =  [x.replace('_GSE176078','') for x in df_umap[df_umap['label']=='GSE176078']['cell']]
dflabel = pd.merge(dflabel,dfl,right_on='cell',left_on='l2',how='left')

print(dfl.columns)
label_index=8
label = dfl.columns[label_index]

df_umap = pd.merge(df_umap,dflabel[['l1',label]],right_on='l1',left_on='cell',how='left')
df_umap[label] = df_umap[label].mask(df_umap[label].isna(), df_umap['label'])
cp = sns.color_palette(cc.glasbey, n_colors=len(df_umap[label].unique()))
p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue=label,s=1,palette=cp)
plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
# t.axes.set_title("UMAP plot of cell-mixture",fontsize=30)
p.set_xlabel("UMAP1",fontsize=20)
p.set_ylabel("UMAP2",fontsize=20)
plt.savefig(spruce.model_id+'_ct_umap_celltype_label.png');plt.close()

##################################################################

# cluster cells using kmeans for cell annotation
from sklearn.cluster import KMeans

dfh = spruce.cell_topic.h.copy()

kmeans = KMeans(n_clusters=50, random_state=0).fit(dfh.iloc[:,1:].to_numpy())

df_umap['cluster'] = kmeans.labels_

cp = sns.color_palette(cc.glasbey, n_colors=len(df_umap['cluster'].unique()))
p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='cluster',s=1,palette=cp,legend=False)

# plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
# t.axes.set_title("UMAP plot of cell-mixture",fontsize=30)
p.set_xlabel("UMAP1",fontsize=20)
p.set_ylabel("UMAP2",fontsize=20)
plt.savefig(spruce.model_id+'_ct_h_kmeans.png');plt.close()


##################################################################

# cell annotation based on kmeans cluster and majority 

df_umap['label'] = [x.split('_')[len(x.split('_'))-1] for x in df_umap['cell']]
f='/home/BCCRC.CA/ssubedi/projects/data/GSE176078mix/GSE176078_metadata.csv.gz'

dfl = pd.read_csv(f,compression='gzip')
dfl = dfl.rename(columns={'Unnamed: 0':'cell'})

dflabel = pd.DataFrame()
dflabel['l1'] =  [x for x in df_umap[df_umap['label']=='GSE176078']['cell']]
dflabel['l2'] =  [x.replace('_GSE176078','') for x in df_umap[df_umap['label']=='GSE176078']['cell']]
dflabel = pd.merge(dflabel,dfl,right_on='cell',left_on='l2',how='left')

print(dfl.columns)
label_index=8
label = dfl.columns[label_index]

df_umap = pd.merge(df_umap,dflabel[['l1',label]],right_on='l1',left_on='cell',how='left')
df_umap[label] = df_umap[label].mask(df_umap[label].isna(), df_umap['label'])

df_clust_celltype = df_umap.groupby(['cluster','celltype_major'])['cell'].count().reset_index().sort_values(['cluster','cell'])
df_clust_celltype = df_clust_celltype.drop_duplicates(subset=['cluster'], keep='last').drop(['cell'],axis=1).rename(columns={'celltype_major':'cluster_celltype'})

df_umap = pd.merge(df_umap,df_clust_celltype,on='cluster',how='left')

cp = sns.color_palette(cc.glasbey, n_colors=len(df_umap['cluster_celltype'].unique()))
p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='cluster_celltype',s=1,palette=cp,legend=True)

# plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
# t.axes.set_title("UMAP plot of cell-mixture",fontsize=30)
p.set_xlabel("UMAP1",fontsize=20)
p.set_ylabel("UMAP2",fontsize=20)
plt.savefig(spruce.model_id+'_ct_h_kmeans_celltype.png');plt.close()

df_umap.to_csv(spruce.model_id+'_ct_h_umap_cordinates_kmeans.csv.gz',index=False,compression='gzip')


##################################################################
                 ## histogram of h ##
##################################################################

dfh = spruce.cell_topic.h
df_hsum = dfh.iloc[:,1:].sum(0).reset_index().rename(columns={'index':'topic',0:'sum'})
df_hsum = df_hsum.sort_values('sum',ascending=False)
p = sns.barplot(x='topic',y='sum',data=df_hsum,color='blue')
p.set_xlabel("Topic",fontsize=20)
p.set_ylabel("Sum",fontsize=20)
plt.savefig(spruce.model_id+'_ct_h_sum.png');plt.close()

df_hmax = pd.DataFrame(pd.Series([x.replace('h','') for x in dfh.iloc[:,1:].idxmax(axis=1)]).value_counts()).reset_index().rename(columns={'index':'topic',0:'argmax_count'})
p = sns.barplot(x='topic',y='argmax_count',data=df_hmax,color='blue')
p.set_xlabel("Topic",fontsize=20)
p.set_ylabel("Count(argmax)",fontsize=20)
plt.savefig(spruce.model_id+'_ct_h_argmax.png');plt.close()

##################################################################


from analysis import _topics

top_n = 5
df_top_genes = _topics.topic_top_genes(spruce,top_n)
df_top_genes.to_csv(spruce.model_id+'_ct_beta_weight_top_'+str(top_n)+'_genes.csv.gz',index=False,compression='gzip')


df_h_sample = _topics.sample_cells_with_celltype(spruce)
df_h_sample.to_csv(spruce.model_id+'_ct_h_sample_celltype.csv.gz',index=False,compression='gzip')


# get previously saved umap coordinates and kmeans cell type label

dfh = spruce.cell_topic.h.copy()
df_umap = pd.read_csv(spruce.model_id+'_ct_h_umap_cordinates_kmeans.csv.gz',compression='gzip')

df_h_sample_kmeans = df_umap.groupby('cluster_celltype').sample(n=50, random_state=1)
dfh = pd.merge(dfh,df_h_sample_kmeans[['cell','cluster_celltype']],on='cell',how='inner')
dfh = dfh.rename(columns={'cluster_celltype':'Topic'})
dfh.to_csv(spruce.model_id+'_ct_h_sample_kmeans_celltype.csv.gz',index=False,compression='gzip')
del dfh

df_h_sample_kmeans = df_umap.groupby('cluster').sample(n=50, random_state=1)
dfh = pd.merge(dfh,df_h_sample_kmeans[['cell','cluster']],on='cell',how='inner')
dfh = dfh.rename(columns={'cluster':'Topic'})
dfh.to_csv(spruce.model_id+'_ct_h_sample_kmeans_cluster.csv.gz',index=False,compression='gzip')
del dfh


