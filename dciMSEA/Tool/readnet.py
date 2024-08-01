import pandas as pd
import numpy as np
import xlrd

def readnet(loadname):

    file_location4 = loadname
    data = xlrd.open_workbook(file_location4)
    table = data.sheets()[0]
    nrows = table.nrows
    ncols = table.ncols
    basenet = np.zeros((nrows, ncols))
    for x in range(ncols):
        cols = table.col_values(x)
        cols1 = np.matrix(cols)
        basenet[:, x] = cols1

    return  basenet