import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import numpy as np
import sys
import pandas as pd
import os

yCells = 1
xCells = 100
numSpecs = 1

def readData(fname):
    results = {}
    with open(fname, 'r') as f:
        key = None
        rowCount = 0
        for i, line in enumerate(f.readlines()):
            splitline = line.split()
            if line.startswith("Time"):
                rowCount = 0
                dataArray = np.zeros((yCells*xCells,numSpecs+2))
                key = splitline[1]
                results.update({splitline[1]: dataArray})

            else:
                dataArray = results[key]
                for col, val in enumerate(splitline):
                    dataArray[rowCount,col] = float(val)
                rowCount += 1

    return results


fname = sys.argv[1]
results = readData(fname)
cwd = os.getcwd()

for imageCount, key in enumerate(results.keys()):
    data = results[key]
    dataArray = np.zeros((yCells,xCells))
    for row in xrange(data.shape[0]):
        dataArray[int(data[row,1]), int(data[row,0])] =  data[row,2]
    print "Image: ",imageCount, " out of: ", len(results.keys())-1
    ax = sns.heatmap(dataArray)
    ax.invert_yaxis()
    plt.ylabel("y")
    plt.xlabel("x")
    plt.title('Time: ' + key + ' sec')
    figName = 'spec.png'
    #plt.show()
    plt.savefig(figName, dpi=500)
    plt.close()

