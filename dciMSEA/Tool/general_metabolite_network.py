
import pandas as pd
import os
import pandas as pd
import numpy as np
import networkx as nx


def get_general_metabolite_network(Cmatrixload,Hmatrixload):
    '''

    :param Cmatrixload: Concentration matrix file address for CRC,

                        The first column is the detected metabolite KEGG id, and each subsequent column is the sample

    :param Hmatrixload:  Concentration matrix file address for Control
    :return:  1.N : matrix of Largest connectivity graph of general metabolite network
              2.all_meatbolite_: All metabolites involved in the maximum Connected graph
              3.KEGGid : detected metabolites KEGG id
              4.refdata: Concentration matrix of control group
              5.case: Concentration matrix of CRC group
   '''

    Hmatrix = pd.read_excel(Hmatrixload, header=None, index_col=0)
    Cmatrix = pd.read_excel(Cmatrixload, header=None, index_col=0)

    refdata = Hmatrix.values
    case = Cmatrix.values
    KEGGid = sorted(list(Cmatrix.index))



    # how much pathway hitted

    filepath = "/code/Dataset/path"
    record_d = []
    for i, j, pathlist in os.walk(filepath):
        for p in pathlist:

            pname = filepath + '/' + p
            pathpd = pd.read_excel(pname)
            pm = []  # 每条通路的化合物数(原始网络)

            for c in range(len(pathpd)):
                pm.append(pathpd.loc[c, 'substrate'])
                pm.append(pathpd.loc[c, 'product'])
            pm = list(set(pm))

            f = [k for k in KEGGid if k in pm]
            if len(f) >= 1:
                record_d.append(p)

    # print('step 2 down')
    # print(record_d)
    #   construct network

    record_pairs = []
    record_metabolites = []
    for i in record_d:
        pname_ = filepath + '/' + i
        pathpd_ = pd.read_excel(pname_)

        for c in range(len(pathpd_)):

            a = pathpd_.loc[c, 'substrate']
            b = pathpd_.loc[c, 'product']
            record_metabolites.append(a)
            record_metabolites.append(b)
            c = sorted([a, b])
            if c not in record_pairs:
                record_pairs.append(c)

    record_metabolites = sorted(set(record_metabolites))
    record_pairs = sorted(record_pairs)
    # print(record_metabolites)

    M = np.zeros((len(record_metabolites), len(record_metabolites)))

    for i in record_pairs:
        # print(i)
        aa = record_metabolites.index(i[0])
        bb = record_metabolites.index(i[1])
        M[aa, bb] = 1
        M[bb, aa] = 1

    GMG_graph = nx.from_numpy_array(M)


    largest = max(nx.connected_components(GMG_graph), key=len)
    largest_connected_subgraph = GMG_graph.subgraph(largest)

    all_metabolite_ = [record_metabolites[i] for i in largest_connected_subgraph.nodes()]
    all_metabolite_ = list(set(all_metabolite_))
    N = np.zeros((len(all_metabolite_), len(all_metabolite_)))
    for i in largest_connected_subgraph.edges():
        a = record_metabolites[i[0]]
        b = record_metabolites[i[1]]

        aa = all_metabolite_.index(a)
        bb = all_metabolite_.index(b)
        N[aa, bb] = 1


    return N,all_metabolite_,KEGGid,refdata,case
# # how much pathway hitted
#
# filepath  = "F:\dci-MSEA\Dataset\path"
# record_d = []
# for i, j, pathlist in os.walk(filepath):
#     for p in pathlist:
#
#         pname = filepath + '/' + p
#         pathpd = pd.read_excel(pname)
#         pm = []  # 每条通路的化合物数(原始网络)
#
#         for c in range(len(pathpd)):
#             pm.append(pathpd.loc[c, 'substrate'])
#             pm.append(pathpd.loc[c, 'product'])
#         pm = list(set(pm))
#
#         f = [k for k in KEGGid if k in pm]
#         if len(f)>=1:
#             record_d.append(p)
#
# print('step 2 down')
# print(record_d)
# #   construct network
#
#
# record_pairs = []
# record_metabolites = []
# for i in record_d:
#     pname_ = filepath + '/' + i
#     pathpd_ = pd.read_excel(pname_)
#
#     for c in range(len(pathpd_)):
#
#         a = pathpd_.loc[c, 'substrate']
#         b = pathpd_.loc[c, 'product']
#         record_metabolites.append(a)
#         record_metabolites.append(b)
#         c = sorted([a,b])
#         if c  not in record_pairs:
#             record_pairs.append(c)
#
# record_metabolites = sorted(set(record_metabolites))
# record_pairs = sorted(record_pairs)
# print(record_metabolites)
#
# M = np.zeros((len(record_metabolites),len(record_metabolites)))
#
# for i in record_pairs:
#     print(i)
#     aa = record_metabolites.index(i[0])
#     bb = record_metabolites.index(i[1])
#     M[aa,bb] = 1
#     M[bb,aa] = 1
#
# M = pd.DataFrame(M)
# M.index = record_metabolites
# M.columns = record_metabolites
# print(M)