import networkx as nx
import random


# 6.对重构的原始图进一步得最小生成树,去环处理,输入之前要判断是否做过图形解冻
def removing_cycle(G):
    while True:
        try:
            cycle = nx.find_cycle(G)
            ##print('Cycle found')
            ##print(cycle)
            edge = random.sample(cycle, 1)[0]
            ##print('随机移除的边为:',edge)
            G.remove_edge(edge[0], edge[1])  # 移除环
        except:
            break
    ##print('检测是否有环并去环')
    return G