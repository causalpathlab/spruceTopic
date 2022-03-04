import pandas as pd
import numpy as np

def run_component_analysis(x,y,method,title):

    from sklearn import preprocessing
    import matplotlib.pyplot as plt
    import seaborn as sns
    from mpl_toolkits.mplot3d import Axes3D
    plt.rcParams["figure.figsize"] = [15.50, 10.50]
    plt.rcParams["figure.autolayout"] = True



    model = None
    if method=="tsne":
        from sklearn.manifold import TSNE
        model = TSNE(n_components=3,random_state=0)
    elif method=="pca":
        from sklearn.decomposition import PCA
        model = PCA(n_components=3)
        standardizer = preprocessing.StandardScaler()
        x = standardizer.fit_transform(x) 

    components = model.fit_transform(x)
    df_components = pd.DataFrame(data = components,columns = [method+'1', method+'2', method+'3'])
    df_components['sample'] = y

    # fig = plt.figure(figsize=(12, 9))
    # ax = Axes3D(fig)
    # for grp_name, grp_idx in df_components.groupby('sample').groups.items():
    #     x = df_components.iloc[grp_idx][method+"1"]
    #     y = df_components.iloc[grp_idx][method+"2"]
    #     z = df_components.iloc[grp_idx][method+"3"]
    #     ax.scatter(x,y,z, label=grp_name) 
    # ax.legend()
    # ax.set_xlabel(method+"1")
    # ax.set_ylabel(method+"2")
    # ax.set_zlabel(method+"3")
    # plt.savefig("../output/"+method+"_3d.png");plt.close()
    custpal = ['grey','dodgerblue','green','red','orange','brown','purple']
    sns.scatterplot(data=df_components, x=method+"1", y=method+"2", hue="sample",palette=custpal)
    plt.xlabel(method+'1', fontsize=20)
    plt.ylabel(method+'2', fontsize=20)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,prop={"size":20})
    plt.savefig("../output/"+method+"_"+title+".png");plt.close()

