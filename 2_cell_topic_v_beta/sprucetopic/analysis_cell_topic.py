import experiment_cell_topic 
import pandas as pd

experiment_home = '/projects/experiments/spruce_topic/2_cell_topic_v_beta/'
spr = experiment_cell_topic.get_experiment(experiment_home)
print(spr.cell_topic.model_id)

##################################################################
# 1  umap figures of cell topic
##################################################################

# 1.1 generate and save umap coordinates based on softmaxed latent dimensions 

import umap
import umap.plot

dfh = spr.cell_topic.h.copy()
df_umap= pd.DataFrame()
df_umap['cell'] = dfh['cell']

umap_2d = umap.UMAP(n_components=2, init='random', random_state=0,min_dist=0.0,metric='cosine')
proj_2d = umap_2d.fit(dfh.iloc[:,1:])

df_umap[['umap1','umap2']] = umap_2d.embedding_[:,[0,1]]
df_umap['cell_topic'] = [x.replace('h','') for x in dfh.iloc[:,1:].idxmax(axis=1)]
df_umap.to_csv(spr.cell_topic.model_id+'_ct_h_umap_cordinates.csv.gz',index=False,compression='gzip')

##################################################################
# 1.2 
# get previously saved umap coordinates and
# generate argmax label plot

import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [15, 10]
plt.rcParams['figure.autolayout'] = True
import colorcet as cc
import seaborn as sns

df_umap = pd.read_csv(spr.cell_topic.model_id+'_ct_h_umap_cordinates.csv.gz',compression='gzip')

# if we need to select topics with >100 cells
# selected_topic = list(df_umap.topic.value_counts().index[0:32])
# df_umap = df_umap[df_umap['cell_topic'].isin(selected_topic)]

# plot umap with argmax topic

cp = sns.color_palette(cc.glasbey_dark, n_colors=len(df_umap['cell_topic'].unique()))
p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='cell_topic',s=5,palette=cp,legend=False)
# plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
p.axes.set_title("UMAP plot of cell-mixture with argmax cell topic assignment",fontsize=30)
p.set_xlabel("UMAP1",fontsize=20)
p.set_ylabel("UMAP2",fontsize=20)
plt.savefig(spr.cell_topic.model_id+'_ct_umap_h_argmax_label.png');plt.close()

##################################################################
# 1.3 
# get previously saved umap coordinates and
# generate umap plot with celltype label from metadata and annotated files

df_umap = pd.read_csv(spr.cell_topic.model_id+'_ct_h_umap_cordinates.csv.gz',compression='gzip')
df_umap['label'] = [x.split('_')[len(x.split('_'))-1] for x in df_umap['cell']]

# metadata
f='/home/BCCRC.CA/ssubedi/projects/data/GSE176078mix/GSE176078_metadata.csv.gz'
dfl = pd.read_csv(f,compression='gzip')
dfl = dfl.rename(columns={'Unnamed: 0':'cell'})

dflabel = pd.DataFrame()
dflabel['l1'] =  [x for x in df_umap[df_umap['label']=='GSE176078']['cell']]
dflabel['l2'] =  [x.replace('_GSE176078','') for x in df_umap[df_umap['label']=='GSE176078']['cell']]
dflabel = pd.merge(dflabel,dfl,right_on='cell',left_on='l2',how='left')

print(dfl.columns)
label_index=7
label = dfl.columns[label_index]

df_umap = pd.merge(df_umap,dflabel[['l1',label]],right_on='l1',left_on='cell',how='left')
df_umap[label] = df_umap[label].mask(df_umap[label].isna(), df_umap['label'])

# annotated file
normal_f='/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/output/scanpy/gse_normal_annotated_chetah.csv.gz'
df_nlabel =  pd.read_csv(normal_f,compression='gzip')
df_nlabel = df_nlabel[~df_nlabel['pruned_labels'].isna()]
df_nlabel['pruned_labels'] = [ x.replace(' ','_')+'_n' for x in df_nlabel['pruned_labels']]

df_umap = pd.merge(df_umap,df_nlabel[['cell','pruned_labels']],on='cell',how='left')
df_umap.loc[df_umap['pruned_labels'].notna(),label] = df_umap[df_umap['pruned_labels'].notna()]['pruned_labels'].values
df_umap =df_umap.drop(columns=['pruned_labels'])
df_umap.loc[df_umap[label]=='GSE164898',label] = 'Epithelial_n'
df_umap.loc[df_umap[label]=='GSE156728-CD8',label] = 'T_CD8_pan'
df_umap.loc[df_umap[label]=='GSE156728-CD4',label] = 'T_CD4_pan'
df_umap[label] = df_umap[label].mask(df_umap[label].isna(), df_umap['label'])

# fix label names for visualization

df_umap[label] = [ x.replace('-','_').replace(' ','_')  for x in df_umap[label]]

df_umap = df_umap.sort_values(label)

xrep={
'B_cell_n':'B_n',
'B_cells_Memory':'B',
'B_cells_Naive':'B',

'CAF_n':'CAFs_n',
'CAFs_MSC_iCAF_like':'CAFs',
'CAFs_myCAF_like':'CAFs',

'Cancer_Basal_SC':'Cancer',
'Cancer_Cycling':'Cancer',
'Cancer_Her2_SC':'Cancer',
'Cancer_LumA_SC':'Cancer',
'Cancer_LumB_SC':'Cancer',

'Endothelial_ACKR1':'Endothelial',
'Endothelial_CXCL12':'Endothelial',
'Endothelial_Lymphatic_LYVE1':'Endothelial',
'Endothelial_RGS5':'Endothelial',
'Endothelial_n':'Endothelial_n',

'Epithelial_n':'Epithelial_n',
'Luminal_Progenitors':'Epithelial',
'Mature_Luminal':'Epithelial',
'Myoepithelial':'Epithelial',

'Cycling_Myeloid':'Myeloid',
'DCs':'DCs',
'Dendritic_n':'DCs_n',
'Macrophage':'Macrophage',
'Macrophage_n':'Macrophage_n',
'Mast_n':'Mast_n',
'Monocyte':'Monocyte',

'Myofibroblast_n':'Myofibroblast_n',

'NKT_cells':'NK',
'NK_cells':'NK',
'NK_n':'NK_n',

'PVL_Differentiated':'PVL',
'PVL_Immature':'PVL',
'Cycling_PVL':'PVL',

'Plasma_n':'Plasma_n',
'Plasmablasts':'Plasma',

'Cycling_T_cells':'T',
'CD4_T_cell_n':'T_n',
'T_CD4_pan':'T_pan',
'CD8_T_cell_n':'T_n',
'T_CD8_pan':'T_pan',
'T_cells_CD4+':'T',
'T_cells_CD8+':'T',
'reg._T_cell_n':'T_n'

}

xnew={
'B_cell_n':'B_n',
'B_cells_Memory':'B_Memory',
'B_cells_Naive':'B_Naive',

'CAF_n':'CAFs_n',
'CAFs_MSC_iCAF_like':'CAFs_MSC_iCAF_like',
'CAFs_myCAF_like':'CAFs_myCAF_like',

'Cancer_Basal_SC':'Cancer_Basal_SC',
'Cancer_Cycling':'Cancer_Cycling',
'Cancer_Her2_SC':'Cancer_Her2_SC',
'Cancer_LumA_SC':'Cancer_LumA_SC',
'Cancer_LumB_SC':'Cancer_LumB_SC',

'Endothelial_ACKR1':'Endothelial_ACKR1',
'Endothelial_CXCL12':'Endothelial_CXCL12',
'Endothelial_Lymphatic_LYVE1':'Endothelial_Lymphatic_LYVE1',
'Endothelial_RGS5':'Endothelial_RGS5',
'Endothelial_n':'Endothelial_n',

'Epithelial_n':'Epithelial_n',
'Luminal_Progenitors':'Epithelial_Luminal',
'Mature_Luminal':'Epithelial_MLuminal',
'Myoepithelial':'Epithelial_Myoepi',

'Cycling_Myeloid':'Myeloid_Cycling',
'DCs':'DCs',
'Dendritic_n':'DCs_n',
'Macrophage':'Macrophage',
'Macrophage_n':'Macrophage_n',
'Mast_n':'Mast_n',
'Monocyte':'Monocyte',

'Myofibroblast_n':'Myofibroblast_n',

'NKT_cells':'NKT',
'NK_cells':'NK',
'NK_n':'NK_n',

'PVL_Differentiated':'PVL_Differentiated',
'PVL_Immature':'PVL_Immature',
'Cycling_PVL':'PVL_Cycling',

'Plasma_n':'Plasma_n',
'Plasmablasts':'Plasma',

'Cycling_T_cells':'T_Cycling',
'CD4_T_cell_n':'T_CD4_n',
'T_CD4_pan':'T_CD4_pan',
'CD8_T_cell_n':'T_CD8_n',
'T_CD8_pan':'T_CD8_pan',
'T_cells_CD4+':'T_CD4',
'T_cells_CD8+':'T_CD8',
'reg._T_cell_n':'T_reg_n'

}

df_umap['celltype'] = [ xrep[x] for x in df_umap[label]]

df_umap = df_umap.sort_values('celltype')

# plot all labels
cp = sns.color_palette(cc.glasbey, n_colors=len(df_umap['celltype'].unique()))
p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='celltype',s=2,palette=cp)
plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
p.axes.set_title("UMAP plot of cell-mixture with cell type label",fontsize=30)
# plt.figtext(.1,.8,"singleR/CHETAH for normal dataset",fontsize=20)
p.set_xlabel("UMAP1",fontsize=20)
p.set_ylabel("UMAP2",fontsize=20)
plt.savefig(spr.cell_topic.model_id+'_ct_umap_celltype_label_all.png');plt.close()


df_umap.to_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_celllabel.csv.gz',index=False,compression='gzip')

##################################################################
# 1.5 
# cluster cells using kmeans 
# and assign celltype to clusters
from sklearn.cluster import KMeans

df_umap = pd.read_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_celllabel.csv.gz')

dfh = spr.cell_topic.h.copy()

kmeans = KMeans(n_clusters=50, random_state=0).fit(dfh.iloc[:,1:].to_numpy())

df_kmeans = dfh[['cell']]
df_kmeans['cluster'] = kmeans.labels_
df_umap = pd.merge(df_umap,df_kmeans,on='cell',how='left')

df_clust_celltype = df_umap.groupby(['cluster',label])['cell'].count().reset_index().sort_values(['cluster','cell'])
df_clust_celltype = df_clust_celltype.drop_duplicates(subset=['cluster'], keep='last').drop(['cell'],axis=1).rename(columns={label:'cluster_celltype'})

df_umap = pd.merge(df_umap,df_clust_celltype,on='cluster',how='left')

##### save 
df_umap.to_csv(spr.cell_topic.model_id+'_ct_h_umap_cordinates_kmeans.csv.gz',index=False,compression='gzip')
df_umap = pd.read_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_kmeans.csv.gz')
######

cp = sns.color_palette(cc.glasbey, n_colors=len(df_umap['cluster'].unique()))
p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='cluster',s=1,palette=cp,legend=False)

# plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
p.axes.set_title("UMAP plot of cell-mixture kmeans label",fontsize=30)
p.set_xlabel("UMAP1",fontsize=20)
p.set_ylabel("UMAP2",fontsize=20)
plt.savefig(spr.cell_topic.model_id+'_ct_h_kmeans_label.png');plt.close()


df_umap = df_umap.sort_values('cluster_celltype')
cp = sns.color_palette(cc.glasbey_dark, n_colors=len(df_umap['cluster_celltype'].unique()))
p = sns.scatterplot(data=df_umap, x='umap1', y='umap2', hue='cluster_celltype',s=2,palette=cp,legend=True)
plt.legend(title='Cell type',title_fontsize=18, fontsize=14,loc='center left', bbox_to_anchor=(1, 0.5))
p.axes.set_title("UMAP plot of cell-mixture kmeans label",fontsize=30)
p.set_xlabel("UMAP1",fontsize=20)
p.set_ylabel("UMAP2",fontsize=20)
plt.savefig(spr.cell_topic.model_id+'_ct_h_kmeans_celltype.png');plt.close()


##################################################################
# 2 data for struct plot in R
##################################################################

dfh = spr.cell_topic.h.copy()
df_umap = pd.read_csv(spr.cell_topic.model_id+'_ct_h_umap_cordinates_kmeans.csv.gz',compression='gzip')

df_h_sample_kmeans = df_umap.groupby('cluster_celltype').sample(n=50, random_state=1)
dfh = pd.merge(dfh,df_h_sample_kmeans[['cell','cluster_celltype']],on='cell',how='inner')
dfh = dfh.rename(columns={'cluster_celltype':'Topic'})
dfh.to_csv(spr.cell_topic.model_id+'_ct_h_sample_kmeans_celltype.csv.gz',index=False,compression='gzip')
del dfh

# dfh = spr.cell_topic.h.copy()
# df_h_sample_kmeans = df_umap.groupby('cluster').sample(n=50, random_state=1)
# dfh = pd.merge(dfh,df_h_sample_kmeans[['cell','cluster']],on='cell',how='inner')
# dfh = dfh.rename(columns={'cluster':'Topic'})
# dfh.to_csv(spr.cell_topic.model_id+'_ct_h_sample_kmeans_cluster.csv.gz',index=False,compression='gzip')
# del dfh

# dfh = spr.cell_topic.h.copy()
# df_umap = pd.read_csv(spr.cell_topic.model_id+'_ct_h_umap_cordinates_kmeans.csv.gz',compression='gzip')
# df_h_sample = df_umap.groupby(label).sample(n=50, random_state=1)
# dfh = pd.merge(dfh,df_h_sample[['cell',label]],on='cell',how='inner')
# dfh = dfh.rename(columns={label:'Topic'})
# dfh.to_csv(spr.cell_topic.model_id+'_ct_h_sample_celltype.csv.gz',index=False,compression='gzip')
# del dfh

# dfh = spr.cell_topic.h.copy()
# df_umap = pd.read_csv(spr.cell_topic.model_id+'_ct_h_umap_cordinates_kmeans.csv.gz',compression='gzip')
# df_umap2 = df_umap[df_umap[label]=='Cancer_Epithelial']
# # df_umap2 = df_umap[df_umap['cluster_celltype']=='Cancer_Epithelial']
# df_h_sample = df_umap2.groupby('cell_topic').sample(frac=0.2, random_state=1)
# dfh = pd.merge(dfh,df_h_sample[['cell','cell_topic']],on='cell',how='inner')
# dfh = dfh.rename(columns={'cell_topic':'Topic'})
# dfh.to_csv(spr.cell_topic.model_id+'_ct_h_sample_cancer_topic.csv.gz',index=False,compression='gzip')
# del dfh

##################################################################
# 3 histogram of h 
##################################################################

dfh = spr.cell_topic.h
df_hsum = dfh.iloc[:,1:].sum(0).reset_index().rename(columns={'index':'cell_topic',0:'sum'})
df_hsum = df_hsum.sort_values('sum',ascending=False)
p = sns.barplot(x='cell_topic',y='sum',data=df_hsum,color='blue')
p.set_xlabel("Topic",fontsize=20)
p.set_ylabel("Sum",fontsize=20)
plt.savefig(spr.model_id+'_ct_h_sum.png');plt.close()

df_hmax = pd.DataFrame(pd.Series([x.replace('h','') for x in dfh.iloc[:,1:].idxmax(axis=1)]).value_counts()).reset_index().rename(columns={'index':'cell_topic',0:'argmax_count'})
p = sns.barplot(x='cell_topic',y='argmax_count',data=df_hmax,color='blue')
p.set_xlabel("Topic",fontsize=20)
p.set_ylabel("Count(argmax)",fontsize=20)
plt.savefig(spr.model_id+'_ct_h_argmax.png');plt.close()

##################################################################
# 4 data for top genes heatmap in R 
##################################################################
 
from analysis import _topics

top_n = 25
df_top_genes = _topics.topic_top_genes(spr,top_n)
df_top_genes.to_csv(spr.cell_topic.model_id+'_ct_beta_weight_top_'+str(top_n)+'_genes.csv.gz',index=False,compression='gzip')


#######
df_umap = pd.read_csv(spr.cell_topic.model_id +'_ct_h_umap_cordinates_kmeans.csv.gz')

df_umap_sum = df_umap.groupby(['cell_topic','cluster_celltype'])['cell'].count().reset_index()
df_umap_sum = df_umap_sum.sort_values(['cell_topic','cell'])
df_umap_sum = df_umap_sum.drop_duplicates(subset=['cell_topic'], keep='last')
df_umap_sum = df_umap_sum[df_umap_sum['cell']>100]
df_umap_sum = df_umap_sum.sort_values('cluster_celltype')


# ignore = [26,14,45,43,16,4,42,48,17,19,50,30,37,11,9,6]
# df_umap_sum = df_umap_sum[~df_umap_sum['cell_topic'].isin(ignore)]

df_umap_sum = df_umap_sum.sort_values('cell',ascending=False)

# df_umap_sum = df_umap_sum[df_umap_sum['cluster_celltype'].str.contains('Cancer')]

top_n = 25
df_top_genes = _topics.topic_top_genes(spr,top_n)
sel_topics = ['k'+str(x) for x in df_umap_sum['cell_topic']]
df_top_genes = df_top_genes[df_top_genes['Topic'].isin(sel_topics)]

df_top_genes = df_top_genes.pivot(index='Topic',columns='Gene',values='Proportion')
df_top_genes[df_top_genes>20] = 20
df_top_genes[df_top_genes<-20] = -20

topic_label = [ x for x in zip(df_umap_sum['cell_topic'],df_umap_sum['cluster_celltype'])]
tld = {}
for x in topic_label:tld[x[0]]= x[1]+'/' + str(x[0])

df_top_genes.index = [ tld[int(x.replace('k',''))] for x in df_top_genes.index]
df_top_genes.to_csv(spr.cell_topic.model_id+'_ct_beta_weight_top_'+str(top_n)+'_genes_selected.csv.gz',index=True,compression='gzip')

top3=[]
for i in range(df_top_genes.shape[0]):top3.append(df_top_genes.iloc[i,:].sort_values(ascending=False).index[:3].values)

df_top3 = pd.DataFrame(top3,columns=['a','b','c'])
df_top3['genes'] =  df_top3['a'] +'/' +df_top3['b']+'/'+df_top3['c']
df_top3.index = df_top_genes.index
df_top3 = df_top3.drop(columns=['a','b','c'])
df_top3.to_csv(spr.cell_topic.model_id+'_ct_beta_weight_top_'+str(top_n)+'_genes_top3.csv.gz',index=True,compression='gzip')
