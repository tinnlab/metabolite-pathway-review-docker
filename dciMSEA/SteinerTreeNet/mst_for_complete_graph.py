from SteinerTreeNet.highlight_tree import highlight_tree
from SteinerTreeNet.metric_closure import metric_closure
from SteinerTreeNet.show_weight_network import show_weight_network
from SteinerTreeNet.itera import itera
import networkx as nx
from networkx.utils import pairwise, not_implemented_for
import matplotlib.pyplot as plt
from itertools import combinations, chain
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
import random

# 3.感兴趣节点的完全连通图 ——> 最小生成树(不具有唯一性)
'''
G:原始图
terminal_nodes:感兴趣的节点

N:更新权重的完全联通图
re:完全连通图得到的最小生成树的边
cost:最小生成树的总权重
constructG:根据最小生成树得到的原始图,可能含有环
T2.nodes():用于后续 次小生成树 的验证
'''



def mst_for_complete_graph(G, N,terminal_nodes):
    constructG = []  # 收集重构图
    re = []
    cost = 0
    mst_edges = nx.minimum_spanning_edges(N, weight='distance', data=True)  # 最小生成树-1
    L12 = []  # 连通图 最小生成树 边
    L3 = list()  # = deepcopy (u,v,{'distance':1,'path':2}) 中的字典
    for u, v, d in mst_edges:
        cost += d['distance']
        re.append((u, v))
        L12.append((u, v))
        L3.append(d)
    # 迭代后,mst_edges已经为空 []
    # edges = chain.from_iterable(pairwise(t['path']) for t in L3) #只有一种可能,输出可以直接映射到图上
    ##print('*********************')
    ##print('L12:', L12)
    ##print('*******************')
    a = 0
    subedges = []
    for t12 in L12:
        #print('*****************************************')
        paths = nx.all_shortest_paths(G, t12[0], t12[1])
        ##print('连通图最小树的边:', t12[0], t12[1])
        allpaths = []  # 针对某条边,输出任意一个映射路径
        for p in paths:
            #print('最短路径长度:',len(p))
            allpaths.append(p)
        samplepath = random.sample(allpaths, 1)[0]  # 随机选取列表中的元素!!!!!!!随机
        a += len(samplepath)
        iteredge = itera(samplepath)
        subedges += iteredge

    T1 = N.edge_subgraph(L12)
    ##print('T1的边数:', T1.number_of_edges())
    ##plt.title('原始最小生成树')
    ##highlight_tree(T1, terminal_nodes)

    T2 = G.edge_subgraph(subedges)
    ##plt.title('原始映射图')
    ##highlight_tree(T2, terminal_nodes)

    constructG.append(T2)

    return N, cost, re, constructG