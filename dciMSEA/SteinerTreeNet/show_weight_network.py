
import networkx as nx
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False

# 2.画权重图

def  show_weight_network(N):

    pos = nx.spring_layout(N)
    nx.draw(N, pos)
    edge_lables = nx.get_edge_attributes(N, 'distance')
    nx.draw_networkx_edge_labels(N, pos, edge_labels=edge_lables)
    nx.draw_networkx_labels(N, pos, alpha=0.5)

    return plt.show()