import numpy as np
from numpy import linalg as LA
import sys

def readData(fname):
    array = np.zeros((2689, 2689))
    with open(fname, 'r') as f:
        for i, line in enumerate(f.readlines()):
            splitline = line.split()
            for col, val in enumerate(splitline):
                array[i,col] = splitline[col]
    return array

transitionMatrix = readData(sys.argv[1])
transitionMatrix = transitionMatrix
print transitionMatrix[0,-1]
u, s, vh = LA.svd(transitionMatrix)
rank = LA.matrix_rank(transitionMatrix)
print s[0], s[rank-1], s[0]/s[rank-1]
