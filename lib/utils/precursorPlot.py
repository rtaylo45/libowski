import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd
import os

def readData(fname):
    results = {}
    with open(fname, 'r') as f:
        key = None
        for i, line in enumerate(f.readlines()):
            splitline = line.split()
            if line.startswith("Time"):
                dataArray = np.zeros((16,7))
                key = splitline[1]
                results.update({splitline[1]: dataArray})

            else:
                dataArray = results[key]
                row = int(splitline[0])
                col = int(splitline[1])
                value = float(splitline[2])
                dataArray[row,col] = value

    return results

fname = sys.argv[1]
results = readData(fname)
xaxis = range(1,17)

cwd = os.getcwd()
imgd = cwd+'/precursorImages'
os.chdir(imgd)

for key in results.keys():
    data = results[key]
   
    plt.grid()
    for spec in xrange(data.shape[1]-1):
        plt.plot(xaxis,data[::,spec],label='C' + str(spec+1))
    plt.xlabel('Reactor height')
    plt.ylabel('Precursor Concentration')
    plt.legend()
    figName = 'precursorst=_'+key+'.png'
    plt.savefig(figName, dpi=800)


