
from Tool import readnet
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix
from Tool.normalizer import row_normalize
import numpy as np
from Tool.getterminal_edge_id import getterminal_edge_id
from Tool.DC_sample import get_DC_sample,mapedges
from Tool.DE_sample import get_DE_sample
from Tool.SteinerNet import SteinerNet

'''

dci-MSEA for individual specific sample if sample in cancer group

'''





# 1.loaddata
from Tool.general_metabolite_network import get_general_metabolite_network


loadC = "/code/Dataset/testCRC.xlsx"
loadH = "/code/Dataset/testControl.xlsx"

GMN_matrix, all_metabolite, detected_metabolites, Control, CRC  = get_general_metabolite_network(loadC, loadH)
GMG_graph = nx.from_numpy_array(GMN_matrix)
GedgesName = [(all_metabolite[ll[0]], all_metabolite[ll[1]]) for ll in GMG_graph.edges()]





# 2. Differential correlations calculation for a specific sample
i = 5  # i represents the i-th sample
S = get_DC_sample(Control, CRC[:, i], 0.05)
DC = mapedges(S, detected_metabolites)

# 3. Differential expressions calculation for a specific sample

DE = get_DE_sample(Control ,CRC[: ,i])


# 4. Differential correlations-informed metabolite network for a specific sample

from Tool.dci_Net import get_dci_Net

W_ = get_dci_Net(GMG_graph, DC, all_metabolite)

# 5. RW-based MSEA

from Tool.RW_based_MSEA import get_MSEA
filepath = "/code/Dataset/path"
t=3
PA = get_MSEA(W_, all_metabolite, detected_metabolites, DE, filepath, t)
print(PA)





