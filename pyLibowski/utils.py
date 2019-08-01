import numpy as np

def readInputMatrix(fname):

    f = open(fname, 'r')
    inpmatrix = f.readlines()
    numcolumns, numrows = len(inpmatrix), len(inpmatrix)
    retmatrix = np.zeros((numrows,numcolumns))

    for rowindex, row in enumerate(inpmatrix):
        splitrow = row.split()
        for columnindex, matrixval in enumerate(splitrow):
            retmatrix[rowindex, columnindex] = float(matrixval)

    return retmatrix
