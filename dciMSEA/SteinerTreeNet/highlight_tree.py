
import networkx as nx
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False

# 7.高亮显示特定点
def highlight_tree(G,terminalnodes):
    allnodes = G.nodes()
    colormap = []
    pos = nx.spring_layout(G)
    for i in allnodes:
        if i in terminalnodes:
            colormap.append('red')
        else:
            colormap.append('blue')
    nx.draw(G,pos,node_color = colormap,with_labels = True,)
    edge_lables = nx.get_edge_attributes(G, 'distance')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_lables)
    nx.draw_networkx_labels(G, pos, alpha=0.5)
    return plt.show()