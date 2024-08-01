'''

输入：terminal_edges,GG
    * terminal_edges:差异边集合 [(),(),...,()]
    * GG:KEGG代谢反应网络的最大联通图

输出：
    * D:  type dict ,将差异边和对应的随机路径整合起来  {差异边1：最短路径1，差异边2：最短路径2,...}
    * D_diff_nodes:  (差异边对应的)随机路径上的节点集合
    * non_D_diff_nodes:  GG.nodes()中不属于pathnodes的节点集合
    * path:  (已经迭代完的路径)用于子图提取的边映射

'''
from SteinerTreeNet.itera import itera
import random
import networkx as nx
import matplotlib.pyplot as plt
def  st_1(terminal_edges,GG):

    pathnodes = []  # 用于删减节点
    path = [] #具体边  映射
    D_diff = []  # 具体路径
    for e in terminal_edges:
        p = nx.all_shortest_paths(GG, e[0], e[1])  # 随机性  迭代完就没
        pp = list(p)
        random_pp = random.sample(pp, 1)[0]
        D_diff.append(itera(random_pp))
        pathnodes += random_pp
        path += itera(random_pp)

    D = dict(zip(terminal_edges, D_diff))
    pathnodes = list(set(pathnodes))

    allnodes = GG.nodes()
    nonpathnodes = [i for i in allnodes if i not in pathnodes]

    #print('D:差异边的具体路径 ', D)
    #print('path:用于子图提取时候的边映射', path)
    #print('pathnodes:差异路径上的节点集合', pathnodes)

    return D,path,pathnodes,nonpathnodes