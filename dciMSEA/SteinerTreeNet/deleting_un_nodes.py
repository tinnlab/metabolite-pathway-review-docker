

import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations, chain
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False

'''
当叶子节点是斯坦纳点的时候,将其删除
输入:
    *removing_cycle_G:  去环以后的图
    *terminal_nodes:    所需要的斯坦纳点

输出:
    *removing_cycle_G:  不含冗余的图
'''

def deleting_un_nodes(removing_cycle_G, terminal_nodes):
    nodes = removing_cycle_G.nodes()  # 获取所有节点
    n = [i for i in nodes if removing_cycle_G.degree(i) == 1]  # 找出度为1的节点
    rn = [i for i in n if i not in terminal_nodes]  # 找出非terminal_nodes的节点
    ##print('非terminal_nodes的节点:',rn)
    for j in rn:
        removing_cycle_G.remove_node(j)  # 删除

    return removing_cycle_G