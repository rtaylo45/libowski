import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd
import seaborn as sns

def readData(fname):
    results = {}
    with open(fname, 'r') as f:
        for i, line in enumerate(f.readlines()):
            splitline = line.split()
            if line.startswith("Size"):
                n = splitline[1]
                header = str(n)
                results[header] = []
            if line.startswith("Time"):
                results[header].append(float(splitline[1]))
    df = pd.DataFrame(data=results)
    return df

def builddf(frames):
    results = {}
    #headers = ['1', '2', '3', '4', '5', '6', '7', '8']
    headers = ['1', '2', '3', '4']
    for i, frame in enumerate(frames):
        header = headers[i]
        results[header] = frame
    df = pd.DataFrame(data=results)
    return df

def builddfSpeedUp(frames):
    results = {}
    #headers = ['1', '2', '3', '4', '5', '6', '7', '8']
    headers = ['1', '2', '3', '4']
    for i, frame in enumerate(frames):
        header = headers[i]
        results[header] = frame
    for header in headers:
        results[header] = results[header]/results['1']
    df = pd.DataFrame(data=results)
    return df

df1 = readData(sys.argv[1])
df2 = readData(sys.argv[2])
df3 = readData(sys.argv[3])
df4 = readData(sys.argv[4])
#df5 = readData(sys.argv[5])
#df6 = readData(sys.argv[6])
#df7 = readData(sys.argv[7])
#df8 = readData(sys.argv[8])
#dflist = [df1, df2, df3, df4, df5, df6, df7, df8]
dflist = [df1, df2, df3, df4]


for key in df1.keys():
    n = key
    frames = []            
    for df in dflist:      
        frames.append(df[key])
    plotDF = builddf(frames)
                           
                           
    ax = sns.catplot(data = plotDF, kind='violin')
    plt.ylabel("Run Time (sec)")
    plt.xlabel("Number of Processors")
    plt.title("Matirx Size = "+key +'x'+key)
    plt.tight_layout()
    plt.savefig("Size"+key+".png", dpi=500)
    plt.close()

    plotdata = []
    #headers = ['1', '2', '3', '4', '5', '6', '7', '8']
    headers = ['1', '2', '3', '4']
    #for header in headers:
    #    plotdata.append(plotDF['1']/plotDF[header])
    plotDf = (plotDF['1']/plotDF)

    plotDF = builddfSpeedUp(frames)
    #plt.plot(headers, plotdata)
    ax = sns.catplot(data = plotDF, kind='violin')
    plt.ylabel('Speep up')
    plt.xlabel('Number of Processors')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.title("Matirx Size = "+key +'x'+key)
    plt.savefig("SpeedUpSize"+key+".png", dpi=500)
    plt.close()














