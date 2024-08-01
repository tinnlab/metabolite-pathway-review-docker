
import networkx as nx
import random


'''

将'm0'和'm1' 对应的两个最短路径连接起来
需要将两个最短路径之间的最短距离边找出来，任选一个

输入:
    *g: 斯坦纳树
    *G: 最原始的图
    *name:  M_diff
    *D2:    {'m0': [9, 13], 'm1': [1, 3], 'm2': [16, 17, 15]}

输出:
     
'''
def    map_to_originalGraph(g,G,name,D2):

    e = g.edges()
    path = []
    for ee in e:

        if ee[0] in name:
            a = D2[ee[0]]
        else:
            a = [ee[0]]
        if ee[1] in name:
            b = D2[ee[1]]
        else:
            b = [ee[1]]

        #print('a:',a)
        #print('b:',b)
        t = []
        for aa in a:
            for bb in b:
                p = nx.all_shortest_paths(G, aa, bb)
                pp = list(p)[0]
                cost = len(pp)-1
                if cost == 1:
                    t.append(tuple(pp))
        #print('t:',t)
        #print('random.sample(t,1):',random.sample(t,1))
        path.append(random.sample(t,1)[0])

    #print('path:',path)

    return path