
import random
import networkx as nx
import matplotlib.pyplot as plt
'''
找出某条路径的一阶邻居，(该路径作为一个整体)

输入:
    *d: 路径节点
    *G: 最大连通图

输出:
    *neighbors:邻居节点

'''

def  neighbor_list(d,G):
    neigh = []
    GG = nx.Graph(G)  #在copy图上进行操作
    for j in d:
        neigh += list(GG.neighbors(j))
        neigh = list(set(neigh))
    neighbors = [i for i in neigh if i not in d]

    return neighbors