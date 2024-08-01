
import numpy as np

def get_DE_sample(refData, sample):   #75个代谢物z-score归一化
    Z = []
    for i in range(len(sample)):
        a = sample[i]  # 代谢物a的浓度值

        mean = np.sum(refData[i, :]) / refData.shape[1]

        std = np.std(refData[i, :])

        zscore = abs(a - mean) / std

        Z.append(zscore)
    zall = np.sum(Z)
    norm_Z = [x/zall for x in Z]

    return norm_Z# Z,zall
