import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import numpy as np
import sys
import pandas as pd
import os

xCells = 14
yCells = 32
numSpecs = 6
maxs = [0.00411086, 0.0187643, 0.0115526, 0.0152787, 0.00151854, 0.000133036]

def readData(fname):
    results = {}
    with open(fname, 'r') as f:
        key = None
        rowCount = 0
        for i, line in enumerate(f.readlines()):
            splitline = line.split()
            if line.startswith("Time"):
                rowCount = 0
                dataArray = np.zeros((xCells*yCells,numSpecs+2))
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
imgd = cwd+'/precursorMultiChanImages'
os.chdir(imgd)

"""
for specID in xrange(1,2):
    groupDir = imgd+'/Group'+str(specID)
    os.chdir(groupDir)
    for imageCount, key in enumerate(results.keys()):
        data = results[key]
        dataArray = np.zeros((16,14))
        for row in xrange(data.shape[0]):
            dataArray[int(data[row,1]), int(data[row,0])] =  data[row,specID+1]
        print "Image: ",imageCount, " out of: ", len(results.keys())-1
        ax = sns.heatmap(dataArray)
        ax.invert_yaxis()
        plt.ylabel("Core Axial Level")
        plt.xlabel("Core Radial Level")
        plt.title('Time: ' + key + ' sec')
        figName = 'precursorst=_'+key+'.png'
        plt.savefig(figName, dpi=500)
        plt.close()
"""

for imageCount, key in enumerate(results.keys()):
    data = results[key]
    fig, axn = plt.subplots(2,3, sharex=True, sharey=True, figsize=(16,9))
    for axCount, specID in enumerate(xrange(1,7)):
        ax = axn.flat[axCount]
        dataArray = np.zeros((yCells,xCells))
        for row in xrange(data.shape[0]):
            dataArray[int(data[row,1]), int(data[row,0])] =  data[row,specID+1]
        print "Image: ",imageCount, " out of: ", len(results.keys())-1
        h = sns.heatmap(dataArray, ax=ax, yticklabels=2, xticklabels=2, vmin = 0.0, vmax = maxs[specID-1])
        h.invert_yaxis()
        ax.set_title("Group "+str(specID))
        fig.suptitle('Time: ' + key + ' sec')
    fig.text(0.5, 0.04, 'Core Radial Level', ha='center')
    fig.text(0.04, 0.5, 'Core Axial Level', va='center', rotation='vertical')
    figName = 'precursorst=_'+key+'.png'
    plt.savefig(figName, dpi=200)
    plt.close()
#plt.show()
