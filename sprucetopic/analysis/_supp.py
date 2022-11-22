import numpy as np
import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [10, 5]
plt.rcParams['figure.autolayout'] = True
import colorcet as cc
import seaborn as sns
import numpy as np

def plot_marker(spr,df,marker_genes):

    fig, ax = plt.subplots(2,3) 
    ax = ax.ravel()

    for i,g in enumerate(marker_genes):
        if g in df.columns:
            print(g)
            val = np.array([x if x<3 else 3.0 for x in df[g]])
            sns.scatterplot(data=df, x='umap1', y='umap2', hue=val,s=.1,palette="viridis",ax=ax[i],legend=False)

            norm = plt.Normalize(val.min(), val.max())
            sm = plt.cm.ScalarMappable(cmap="viridis",norm=norm)
            sm.set_array([])

            # cax = fig.add_axes([ax[i].get_position().x1, ax[i].get_position().y0, 0.01, ax[i].get_position().height])
            fig.colorbar(sm,ax=ax[i])
            ax[i].axis('off')

            ax[i].set_title(g)
    fig.savefig(spr.cell_topic.id+'5_umap_marker_genes_legend.png',dpi=600);plt.close()
