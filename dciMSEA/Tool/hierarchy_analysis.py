
from readnet import readnet
import numpy as np
import pandas as pd
from sklearn.metrics import jaccard_score
import seaborn as sns
import matplotlib

matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster import hierarchy

def hierarchy_analysis(C, method='average'):
    S = np.zeros((C.shape[1], C.shape[1]))

    for i in range(C.shape[1]):
        for j in range(C.shape[1]):
            jac = jaccard_score(C[:, i], C[:, j])
            S[i, j] = jac
    Spd = pd.DataFrame(S)


    mergings = linkage(Spd, method=method)
    fig = plt.figure(figsize=(5, 2.5))
    Z = mergings
    dflt_col = "#808080"  # Unclustered gray
    D_leaf_colors = {"attr_1": dflt_col,

                     "attr_4": "#B061FF",  # Cluster 1 indigo
                     "attr_5": "#B061FF",
                     "attr_2": "#B061FF",
                     "attr_8": "#B061FF",
                     "attr_6": "#B061FF",
                     "attr_7": "#B061FF",

                     "attr_0": "#61ffff",  # Cluster 2 cyan
                     "attr_3": "#61ffff",
                     "attr_9": "#61ffff",
                     }

    dendrogram(mergings,
               color_threshold=1.72,

               leaf_rotation=0,
               leaf_font_size=5,
               above_threshold_color='black',

               )
    hierarchy.set_link_color_palette(["#B061FF", "#61ffff"])
    plt.xlabel('samples')
    plt.show()

    return mergings

def get_clutsermap(C):

    S = np.zeros((C.shape[1], C.shape[1]))

    for i in range(C.shape[1]):
        for j in range(C.shape[1]):
            jac = jaccard_score(C[:, i], C[:, j])
            S[i, j] = jac
    Spd = pd.DataFrame(S)

    a = sns.clustermap(data=Spd, method = 'average',row_cluster=False,vmin=0.2, vmax=0.5, cmap='YlGnBu', tree_kws={
        #                         #'colors':'steelblue',#线色
        'linewidths': 3})

    plt.show()
    return  a.data2d
