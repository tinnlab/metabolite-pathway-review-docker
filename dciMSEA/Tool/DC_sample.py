

import  numpy as np
import scipy.stats as stats




def get_DC_sample(refData, sample, pTh=0.05):
    '''

    :param refData: reference 1547*16  变量数*样本数
    :param perturedData: single sample data 1547*1  单样本网络的样本
    :param pTh: p value threshold
    :return: SSN Adj(0,1)
    '''

    refPcc = np.corrcoef(refData)


    if np.sum(np.isnan(refPcc)):  # if contains nan turn to 0
        refPcc[np.isnan(refPcc)] = 0
        print('refPcc_Nan_%d' % np.sum(np.isnan(refPcc)))


    perturedData = np.append(refData, sample.reshape(-1, 1), axis=1)  # 变成一列
    print('perturedDta.shape:', perturedData.shape)
    perturedPcc = np.corrcoef(perturedData)
    if np.sum(np.isnan(perturedPcc)):  # if contains nan turn to 0
        refPcc[np.isnan(perturedPcc)] = 0
        print('perturedPcc_Nan_%d' % np.sum(np.isnan(perturedPcc)))

    diffPcc = perturedPcc - refPcc


    N = refData.shape[1]
    zScoreArr = (N - 1) * diffPcc / (1 - refPcc * refPcc + 1e-17)
    zScoreArr[np.diag_indices_from(zScoreArr)] = 0  # exclude duijiao

    if np.sum(np.isnan(zScoreArr)):
        print('z-score nan', np.sum(np.isnan(zScoreArr)))
        zScoreArr[np.isnan(zScoreArr)] = 0
    #
    pValueAArr = np.zeros_like(zScoreArr)
    for row in range(len(zScoreArr)):
        p_value = stats.norm.sf(abs(zScoreArr[row])) * 2
        pValueAArr[row] = p_value

    SSN_r = pValueAArr < pTh  # bool

    if np.sum(SSN_r[np.diag_indices_from(SSN_r)]):  # diag should be 0
        print('SSN_%s diag not 0' % id)

    SSN_r[np.diag_indices_from(SSN_r)] = 0  # make sure diag is 0

    S = SSN_r.astype(int)
    S = np.triu(S)
    #X = zScoreArr
    #X = np.triu(X)
    #X += X.T - np.diag(X.diagonal())
    #print('验证是否对称:',X.T == X)
    #print(S)
    res = []
    for m in range(len(S)):
        for mm in range(len(S)):
            if(S[m,mm] == 1):

                res.append((m,mm))

    return res



def mapedges(idlist,id):

    '''

    :param idlist: nonzero index list
    :param id: name of metabolites
    :return: name of metabolites with nozero index
    '''

    name = [(id[i[0]],id[i[1]]) for i in idlist]


    return  name