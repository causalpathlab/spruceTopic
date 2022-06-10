import pandas as pd
dfz = pd.read_pickle('GSEmixscanpy_raw_pipeline_PC_dimensions.pkl')
import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [12.50, 10.50]
plt.rcParams['figure.autolayout'] = True
import umap
import umap.plot
import igraph
import colorcet as cc
import seaborn as sns
umap_2d = umap.UMAP(n_components=2, init='random', random_state=0,min_dist=0.1)
proj_2d = umap_2d.fit(dfz.iloc[:,1:])
dfz['umap1'] = umap_2d.embedding_[:,0]
dfz['umap2'] = umap_2d.embedding_[:,1]
dfz['label'] = [x.split('_')[len(x.split('_'))-1] for x in dfz['index']]

cp = sns.color_palette(cc.glasbey, n_colors=len(dfz['label'].unique()))

sns.scatterplot(data=dfz, x='umap1', y='umap2', hue='label',s=10,palette=cp)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('test_umap_z_label_topic.png',dpi=300)
plt.close()

f='/home/BCCRC.CA/ssubedi/projects/spruce_topic/input/GSE176078/GSE176078_metadata.csv.gz'
dfl = pd.read_csv(f,compression='gzip')
dfl = dfl.rename(columns={'Unnamed: 0':'cell'})

dflabel = pd.DataFrame()
dflabel['l1'] =  [x for x in dfz[dfz['label']=='GSE176078']['index']]
dflabel['l2'] =  [x.replace('_GSE176078','') for x in dfz[dfz['label']=='GSE176078']['index']]
dflabel = pd.merge(dflabel,dfl,right_on='cell',left_on='l2',how='left')

dfz = pd.merge(dfz,dflabel[['l1','celltype_major']],right_on='l1',left_on='index',how='left')
dfz['celltype_major'] = dfz['celltype_major'].mask(dfz['celltype_major'].isna(), dfz['label'])
cp = sns.color_palette(cc.glasbey, n_colors=len(dfz['celltype_major'].unique()))
sns.scatterplot(data=dfz, x='umap1', y='umap2', hue='celltype_major',s=10,palette=cp)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('test_umap_z_label_celltype.png',dpi=300)
plt.close()
