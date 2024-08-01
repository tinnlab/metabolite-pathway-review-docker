

'''

Algorithm for the differential correlation informed metabolite network construction.


'''
from Tool.readnet import readnet
from Tool.getterminal_edge_id import getterminal_edge_id
from Tool.SteinerNet import SteinerNet
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix
from Tool.normalizer import row_normalize
import numpy as np



def  get_dci_Net(GMG_graph,DC,all_name,):

    ''''

    :param GMG_graph:generla metabolite network graph
    :param DC: differential correlation
    :param all_name:all metabolite name list

    '''
    newG = nx.Graph(GMG_graph)
    DC_id = getterminal_edge_id(DC, all_name)

    stNodes, stEdges = SteinerNet(20, DC_id, newG, newG)

    beta = 1.5
    Dci_Net = nx.Graph(GMG_graph)
    for ii, jj in stEdges:
        Dci_Net.add_edge(ii, jj, weight=Dci_Net[ii][jj]['weight'] * beta)  # w = 1.3

    W = csr_matrix(nx.to_numpy_matrix(Dci_Net))  # W = the adjacency matrix of the Dci-Net
    W_ = row_normalize(W)  # W_ = row-normalized matrix of W
    W_ = W_.T

    return W_
