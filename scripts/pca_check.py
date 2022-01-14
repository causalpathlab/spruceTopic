import pandas as pd
import numpy as np

def run_pca(x,y):
    from sklearn.decomposition import PCA
    from sklearn import preprocessing
    import matplotlib.pyplot as plt
    import seaborn as sns
    from mpl_toolkits.mplot3d import Axes3D

    standardizer = preprocessing.StandardScaler()
    x = standardizer.fit_transform(x) 

    pca = PCA(n_components=3)
    principal_components = pca.fit_transform(x)
    df_principal = pd.DataFrame(data = principal_components,columns = ['pc1', 'pc2', 'pc3'])
    df_principal['sample'] = y

    fig = plt.figure(figsize=(12, 9))
    ax = Axes3D(fig)
    for grp_name, grp_idx in df_principal.groupby('sample').groups.items():
        x = df_principal.iloc[grp_idx]["pc1"]
        y = df_principal.iloc[grp_idx]["pc2"]
        z = df_principal.iloc[grp_idx]["pc3"]
        ax.scatter(x,y,z, label=grp_name) 
    ax.legend()
    ax.set_xlabel("pc1")
    ax.set_ylabel("pc2")
    ax.set_zlabel("pc3")
    plt.savefig("../output/pca_3d.png");plt.close()

    sns.scatterplot(data=df_principal, x="pc1", y="pc2", hue="sample")
    plt.savefig("../output/pca_2d.png");plt.close()

