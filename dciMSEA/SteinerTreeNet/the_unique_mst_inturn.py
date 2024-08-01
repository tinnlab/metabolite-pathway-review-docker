from SteinerTreeNet.highlight_tree import highlight_tree
from SteinerTreeNet.metric_closure import metric_closure
from SteinerTreeNet.show_weight_network import show_weight_network
from SteinerTreeNet.itera import itera
import networkx as nx
from itertools import combinations, chain
from networkx.utils import pairwise, not_implemented_for
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
import random
# 4.确定最小生成树的唯一性 ——> 依次删除边
'''
G1: 感兴趣节点的完全连通图
G2: 原始图
edge_for_rm:找到的一个最小生成树的边list
cost:最小生成树的权重
T1:找到的其他最小生成树
T2:根据T1的重构图
'''
def the_unique_mst_inturn(G1,G2, edge_for_rm, cost,terminal_nodes):
    constructG = []
    #图形解冻，依次删除边
    for i in edge_for_rm:
        H = nx.Graph(G1)
        H.remove_edge(i[0], i[1])
        c = 0
        mst_edges = nx.minimum_spanning_edges(H, weight='distance', data=True)
        L12 = []
        L3 = list()
        for u, v, d in mst_edges:
            c += d['distance']
            L12.append((u,v))
            L3.append(d)
        ##print('*********************')
        ##print('L12:',L12)
        ##print('*******************')
        if c == cost:
            T1 = G1.edge_subgraph(L12)
            edges = chain.from_iterable(pairwise(t['path']) for t in L3)

            subedges = []
            for t12 in L12:
                paths = nx.all_shortest_paths(G2, t12[0], t12[1])
                allpaths = []  # 针对某条边,输出任意一个映射路径
                for p in paths:
                    allpaths.append(p)
                samplepath = random.sample(allpaths, 1)[0]  # 随机选取列表中的元素!!!!!!!随机
                iteredge = itera(samplepath)
                subedges += iteredge

            T2 = G2.edge_subgraph(subedges)

            ##print('find:'+'依次删除边',i[0],i[1])
            ##plt.title('删除边得到的最小生成树'+'  '+'{}'.format([i[0],i[1]]))
            ##highlight_tree(T1,terminal_nodes)
            ##plt.title( '原始图  ' + '{}'.format([i[0], i[1]]))
            ##highlight_tree(T2,terminal_nodes)
            constructG.append(T2)
    return constructG
