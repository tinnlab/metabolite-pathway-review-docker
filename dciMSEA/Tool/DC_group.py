import math
import numpy as np
from scipy.linalg import expm,logm


def get_DC_group(Y,X,nY,nX,KEGGid):

    PCCY = np.corrcoef(Y)
    PCCX = np.corrcoef(X)
    f1 = 0.5*(np.log((1+PCCY)/(1-PCCY+1e-17)))
    f2 = 0.5*(np.log((1+PCCX)/(1-PCCX+1e-17)))

    f3 = math.sqrt(1/(nX-3) + 1/(nY-3))

    DifferZ = (f1-f2)/(f3+1e-17)
    DZ = np.triu(DifferZ)
    DZ = abs(DZ)

    t = []
    for i in range(len(DZ)):
        for j in range(len(DZ)):
            if (DZ[i, j] > 1.96):
                t.append((i, j))


    t1 = [tt for tt in t if tt[0] != tt[1]]
    DC = mapedges(t1, KEGGid)

    return DC



def mapedges(idlist,id):

    '''

    :param idlist: nonzero index list
    :param id: name of metabolites
    :return: name of metabolites with nozero index
    '''

    name = [(id[i[0]],id[i[1]]) for i in idlist]


    return  name