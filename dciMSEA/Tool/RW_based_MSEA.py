
from Tool.normalizer import row_normalize
import numpy as np
from Tool.getterminal_edge_id import getterminal_edge_id
from Tool.getPA import getPA

def get_MSEA(W_,all_name,KEGGid,DE,filepath,t):

    '''

    :param W_: W_ = row-normalized matrix of W
    :param all_name: all metabolite name list
    :param KEGGid:
    :param DE: Differential expression
    :param filepath: directory of pathways list
    :param t: iteration step
    :return: PA: Pathway activity for each pathway
    '''
    NA = np.zeros(W_.shape[1])               # Node activity
    seeds = [ii for ii, xx in enumerate(all_name) if xx in KEGGid]
    NA[seeds] = DE


    for rtt in range(t):
        W_ = np.dot(W_, W_)
    W_ = row_normalize(W_)
    W_ = W_.T

    NA = W_.dot(NA)


    filepath = "/code/Dataset/path"
    # from getPA import getPA
    PA = getPA(filepath,all_name,NA)         # Pathway activity for each pathway


    return PA