import numpy as np
import math


def get_DE_group(Y,X):

    '''

    :param Y: control ,meatbolite j
    :param X: case,metabolite j
    :return:
    '''
    ES = []
    for i in range(Y.shape[0]):

        uY = np.mean(Y[i,:])
        uX = np.mean(X[i,:])

        nY = Y.shape[1]
        nX = X.shape[1]
        f1 = (nX-1)*np.var(X[i,:])+(nY-1)*np.var(Y[i,:]) #  算数平均值   方差
        f2 = nX + nY - 2
        pooled = math.sqrt((f1/f2)*(1/nX + 1/nY))

        k = 1-3/(4*(nX+nY)-9)
        es = abs(((uY - uX)*k)/pooled)

        ES.append(es)



    #esall = np.sum(ES)
    #res = [j/esall for j in ES]


    return   ES #res,np.sum(ES)