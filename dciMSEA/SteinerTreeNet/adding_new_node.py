

import networkx as nx

'''
删除差异边路径的节点，并将整个差异边当作一个新顶点；

输入:
    *k:需要删除的节点
    *neighbor:用于新节点的连接
    *G:已经删除过节点的边，后面在它的基础上继续传递下去；

同时还要检查邻居节点有没有在其他的差异边中

'''
def adding_new_node(k,neighbor,G,a):

    GG = nx.Graph(G)
    for j in k:
        GG.remove_node(j)
    for i in neighbor:
        GG.add_edge(a, i)

    return GG