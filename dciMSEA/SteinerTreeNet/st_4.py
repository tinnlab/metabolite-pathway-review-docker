from SteinerTreeNet.adding_new_node import adding_new_node
from SteinerTreeNet.neighbor_list import neighbor_list
from SteinerTreeNet.adjacement import adjacement,construct_g
import networkx as nx
import matplotlib.pyplot as plt


'''
输入:
    *D2:   zip(name,k) {'m0': [9, 12], 'm1': [4, 8, 11, 7]}
    *k:    合并的分量 [[9, 12], [4, 8, 11, 7]]
    *GG:    原始图
    *name:  新命名的节点 {'m0','m1',...}
g:
'''

def st_4(k,name,GG,D2):
    g = nx.Graph(GG)
    for i in range(len(k)):
        neighbor = neighbor_list(k[i], GG)
        g = adding_new_node(k[i], neighbor, g, name[i])

    gg = construct_g(name, D2, g)

    #print('gg.nodes():', gg.nodes())

    return gg