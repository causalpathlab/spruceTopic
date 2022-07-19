import numpy as np
import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [10, 5]
plt.rcParams['figure.autolayout'] = True
import colorcet as cc
import seaborn as sns
import numpy as np

def plot_marker(spr,df,marker_genes):

    fig, ax = plt.subplots(2,4) 
    ax = ax.ravel()

    for i,g in enumerate(marker_genes):
        if g in df.columns:
            print(g)
            sns.scatterplot(data=df, x='umap1', y='umap2', hue=df[g].values,s=0.1,palette="viridis",legend=False,ax=ax[i])
            # ax[i].set_xlabel("UMAP1",fontsize=20)
            # ax[i].set_ylabel("UMAP2",fontsize=20)
            ax[i].set_title(g)
    fig.savefig(spr.cell_topic.model_id+'_umap_marker_genes.png');plt.close()
