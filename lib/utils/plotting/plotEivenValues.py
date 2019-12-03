import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import sys
import pandas as pd

def readData(fname):
    realParts = []
    imagParts = []
    with open(fname, 'r') as f:
        for i, line in enumerate(f.readlines()):
            splitline = line.split()
            realParts.append(float(splitline[0]))
            imagParts.append(float(splitline[1]))
    return [realParts, imagParts]

cNom = readData(sys.argv[1])
c50PercentNom = readData(sys.argv[2])
c150PercentNom = readData(sys.argv[3])
plotList = [cNom, c50PercentNom]
ax = plt.subplot(111)
plt.scatter(cNom[0], cNom[1], marker='x', c='k', label='Nominal\nFlow Rate')
plt.scatter(c50PercentNom[0], c50PercentNom[1], marker='x', c='r', label='50% Nominal\nFlow Rate')
plt.scatter(c150PercentNom[0], c150PercentNom[1], marker='x', c='b', label='150% Nominal\nFlow Rate')
plt.ylabel("Imaginary")
plt.xlabel("Real")
#plt.ylim(-1.e-3, 1.e-3)
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
#plt.yscale('log')
plt.grid()
plt.tight_layout()
#plt.show()
plt.savefig("eigenvalues.png")

