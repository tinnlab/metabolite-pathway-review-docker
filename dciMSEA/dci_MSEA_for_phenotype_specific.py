
from Tool import readnet
from Tool.getterminal_edge_id import getterminal_edge_id
from Tool.normalizer import row_normalize
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd






# 1.loaddata
from Tool.general_metabolite_network import get_general_metabolite_network


loadC = "/code/Dataset/testCRC.xlsx"
loadH = "/code/Dataset/testControl.xlsx"

GMN_matrix, all_metabolite, detected_metabolites, Control, CRC  = get_general_metabolite_network(loadC, loadH)
GMG_graph = nx.from_numpy_array(GMN_matrix)   # graph form of largest connectivity of generla metabolite network
GedgesName = [(all_metabolite[ll[0]], all_metabolite[ll[1]]) for ll in GMG_graph.edges()]   # metabolite name corresponding to edges


print('down')

# 2. Differential correlations calculation for the POI group

from Tool.DC_group import get_DC_group
DC = get_DC_group(Control, CRC, Control.shape[1], CRC.shape[1], detected_metabolites)

print('down')
# 3. Differential expressions calculation for POI group

from Tool.DE_group import get_DE_group
DE = get_DE_group(Control, CRC)

print('down')
# 4. Differential correlations-informed metabolite network

from Tool.dci_Net import get_dci_Net

W_ = get_dci_Net(GMG_graph, DC, all_metabolite)  # W_ = row-normalized matrix of W,  W = the adjacency matrix of the Dci-Net

print('down')
# 5. RW-based MSEA

from Tool.RW_based_MSEA import get_MSEA
filepath = "/code/Dataset/path"
t = 3
PA = get_MSEA(W_, all_metabolite, detected_metabolites, DE, filepath,t) #Pathway activity for each pathway
print('down')
print(PA)