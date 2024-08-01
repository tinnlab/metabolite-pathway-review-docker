



import os
import pandas as pd

def getPA(filepath,larMRN_name,x):
    allpscore = {}
    allscore = []

    for i, j, pathlist in os.walk(filepath):
        for p in pathlist:

            pscore = 0

            pname = filepath + '/' + p
            pathpd = pd.read_excel(pname)
            each_path_name = []  # 每条通路的化合物数(原始网络)


            for c in range(len(pathpd)):
                each_path_name.append(pathpd.loc[c, 'substrate'])
                each_path_name.append(pathpd.loc[c, 'product'])


            each_path_name_inLargest = list(set([n for n in each_path_name if n in larMRN_name] )) # 每条通路的化合物数(最大图)

            epni_index = [y for y,q in enumerate(larMRN_name)if q in each_path_name_inLargest]

            for k in epni_index:
                pscore += x[k]
            allpscore[p] = pscore
            allscore.append(pscore)
            #print('通路分数:',p,pscore)

    return allscore