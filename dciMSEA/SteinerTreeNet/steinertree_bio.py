import networkx as nx
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False

from metric_closure import metric_closure
from mst_for_complete_graph import mst_for_complete_graph
from the_unique_mst_inturn import the_unique_mst_inturn
from the_unique_mst_overall import the_unique_mst_overall

from st_1 import st_1
from st_2 import st_2
from st_3 import st_3
from st_4 import st_4

from removing_cycle import removing_cycle
from deleting_un_nodes import deleting_un_nodes
from map_to_originalGraph import map_to_originalGraph

from tqdm import tqdm

'''
输入:
    *terminal_edges:差异边
    *GG:最大联通子图
    *G:原始图

输出:
    *gIV:斯坦纳树映射后的图

'''

def steinertree_bio(number,terminal_edges,GG,G):
    stN = []
    stE = []


    for i in range(number):
        print('steiner tree num:',i+1)

        D, path, D_diff_nodes, non_D_diff_nodes = st_1(terminal_edges, GG)
        #print('path:', path)
        # k,name,D2 ,= st_2(path)
        setD_diff, M_diff, D2 = st_2(path)
        mapping = st_3(setD_diff, non_D_diff_nodes)
        gg = st_4(setD_diff, M_diff, GG, D2)

        terminal_nodes = M_diff
        M = metric_closure(gg, weight='weight')
        N = M.subgraph(terminal_nodes)
        NN, cost, re, constructG1 = mst_for_complete_graph(gg, N, terminal_nodes)
        constructG2 = the_unique_mst_inturn(NN, gg, re, cost, terminal_nodes)
        constructG3 = the_unique_mst_overall(NN, gg, re, cost, terminal_nodes)
        Glist = constructG1 + constructG2 + constructG3

        for i in Glist:  # 将所有的图依次进行去环
            gI = nx.Graph(i)  # 图形解冻
            gII = removing_cycle(gI)  # 去环
            gIII = deleting_un_nodes(gII, terminal_nodes)  # 删除非必要节点   最终的图
            q = map_to_originalGraph(gIII, G, M_diff, D2)
            subedges = path + q
            #print('q:', q)
            gIV = G.edge_subgraph(subedges)

            stN += gIV.nodes()

            stE += gIV.edges()
    stNodes = (GG.subgraph(stN)).nodes()
    stEdges = (GG.edge_subgraph(stE)).edges()

    return  stNodes, stEdges,GG.edge_subgraph(stE)


def steinertree_bio_u(number,terminal_edges,GG,G):   # 改为并集
    stN = []
    stE = []

    pbar = tqdm(total=number)
    old_elen = 0

    resNodes = []
    resEdges = []
    for i in range(number):
        #print('------斯坦纳树迭代次数-----:',i+1)  20220515 改动
        D, path, D_diff_nodes, non_D_diff_nodes = st_1(terminal_edges, GG)
        #print('path:', path)
        # k,name,D2 ,= st_2(path)
        setD_diff, M_diff, D2 = st_2(path)
        mapping = st_3(setD_diff, non_D_diff_nodes)
        gg = st_4(setD_diff, M_diff, GG, D2)
        #nx.draw(gg, with_labels=True)
        #plt.show()
        terminal_nodes = M_diff
        M = metric_closure(gg, weight='weight')
        N = M.subgraph(terminal_nodes)
        NN, cost, re, constructG1 = mst_for_complete_graph(gg, N, terminal_nodes)
        constructG2 = the_unique_mst_inturn(NN, gg, re, cost, terminal_nodes)
        constructG3 = the_unique_mst_overall(NN, gg, re, cost, terminal_nodes)
        Glist = constructG1 + constructG2 + constructG3

        for ii in Glist:  # 将所有的图依次进行去环
            gI = nx.Graph(ii)  # 图形解冻
            gII = removing_cycle(gI)  # 去环
            gIII = deleting_un_nodes(gII, terminal_nodes)  # 删除非必要节点   最终的图
            q = map_to_originalGraph(gIII, G, M_diff, D2)
            subedges = path + q
            #print('q:', q)
            gIV = G.edge_subgraph(subedges)

            stN += gIV.nodes()

            stE+= gIV.edges()

        stNodes = (GG.subgraph(stN)).nodes()
        stEdges = (GG.edge_subgraph(stE)).edges()


        if(len(stEdges) == old_elen):

            resNodes += stNodes
            resEdges += stEdges
            #print(resEdges)
            pbar.set_description("The  steiner tree iteration has converged at %d-iter" % (i + 1))

            break
        old_elen = len(stEdges)
        pbar.update(1)
    pbar.close()


        # stEdges = (GG.edge_subgraph(stE)).edges(data = True)  20220615
    #stNodes = (GG.subgraph(stN)).nodes()  20220615
    #stEdges = (GG.edge_subgraph(stE)).edges(data = True)  20220615


    return  resNodes,resEdges
