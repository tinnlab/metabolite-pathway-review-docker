

import random
import networkx as nx
import matplotlib.pyplot as plt

'''
判断两个差异路径边是否相邻
方法:邻居节点 vs 差异边路径节点
b:所有邻居组的数组集合
a:新的节点命名字典
g:最后的图(已经有新命名节点的图)
'''
def  adjacement(a,b,D2,g):

    neigh = list(g.neighbors(a))
    #print('adjacement中的neigh:',neigh)
    c = [i for i in neigh if i in D2[b]]
    #print('adjacement中的c:',c)
    d1 = [] #删除
    d2 = [] #连接
    if c!=[]:
        d2 = (a,b)
        for j in c:
            d1 += (a,j)

    return d1,d2

'''
a:命名列表的长度 len(names)
D2:命名节点 vs 差异边节点
g:最后的图(已经有新命名节点的图) 

'''


def construct_g (name,D2,g):  #检测到相邻的边 删除 /连接
    for i in range(0, len(name)):
        for j in range(i + 1, len(name)):
            d1,d2 = adjacement(name[j], name[i], D2, g)
            if (d1 != []) & (d2 != []):

                g.remove_edge(d1[0],d1[1])
                g.add_edge(d2[0],d2[1])
                g.remove_node(d1[1])

    return g