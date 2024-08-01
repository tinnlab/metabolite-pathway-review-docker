import random
import networkx as nx
import matplotlib.pyplot as plt


'''
输入: path ----根据路径迭代的边
输出:
    * k:    每个分量的节点(或 最短路径的节点集)
    * name: 每个分量当做一个整体，命名为一个新节点 {m0,m1,m2,....,}
    * D2:   dict type ,将新命名的节点和分量整合起来，{m0:k[0],m1:k[1],...,}

'''

def st_2(path):

    newG = nx.Graph()  # 创建新图
    newG.add_edges_from(path)

    #newG_conn = list(nx.connected_component_subgraphs(newG))
    newG_conn = list(newG.subgraph(c) for c in nx.connected_components(newG))

    k = []
    for i in newG_conn:
        k.append(list(i.nodes()))

    name = []
    for i in range(len(newG_conn)):
        name.append('m' + str(i))

    D2 = dict(zip(name, k))

    #print('st_2:将每一条最短路径看成新的节点:', name)  # 如果节点只有一个，说明最短路径重合在一起


    return k,name,D2,