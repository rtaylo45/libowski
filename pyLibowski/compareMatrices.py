from __future__ import print_function
import numpy as np
import scipy.linalg as LA
import sys
import os
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

def readInputVector(fname):

    f = open(fname, 'r')
    inpvect = f.readlines()
    numrows = len(inpvect)
    retvect = np.zeros((numrows))

    for rowindex, vectVal in enumerate(inpvect):
        retvect[rowindex] = float(vectVal)

    return retvect

def calcDNNZ(matrix):
    count = 0
    nRows = matrix.shape[0]

    for rowIndex in xrange(nRows):
        diagVal = matrix[rowIndex,rowIndex]
        if diagVal > 0.0:
            count += 1
    return count

def calcONNZ(matrix):
    count = 0
    nRows = matrix.shape[0]
    nColumns = matrix.shape[1]

    for rowIndex in xrange(nRows):
        for columnIndex in xrange(nColumns):
            if rowIndex == columnIndex:
                pass
            else:
                val = matrix[rowIndex,columnIndex]
                if val > 0.0:
                    count += 1
    return count
cwd = os.getcwd()

AMatrix = sys.argv[1]
bVector = sys.argv[2]
solVector = sys.argv[3]

matrixBaseCase = readInputMatrix(AMatrix)

bVectorBaseCase = readInputVector(bVector)

solVectorBaseCase = readInputVector(solVector)

scipySolution = LA.solve(matrixBaseCase, bVectorBaseCase)
print (scipySolution)

print (calcDNNZ(matrixBaseCase))
print (calcONNZ(matrixBaseCase))

k = 0
"""
for i in xrange(0,769, 8):
    c1Row = i
    c2Row = i + 1
    c3Row = i + 2
    c4Row = i + 3
    c5Row = i + 4
    c6Row = i + 5
    xenonRow = i + 6
    iodineRow = i + 7
    iodineSolValBaseCase = solVectorBaseCase[iodineRow]
    #print (np.isclose(nativeSolutionBaseCase[xenonRow],nativeSoltuionPre[xenonRow]))

    # Test the A matrix
    xenonRowPrecursorSourceSetToZero = matrixPrecursorSourceSetToZero[xenonRow,:]
    xenonRowBaseCase = matrixBaseCase[xenonRow,:]
    # checks the row 
    for columnIndex, val in enumerate(xenonRowPrecursorSourceSetToZero):
        testVal = xenonRowBaseCase[columnIndex]
        if testVal != val:
            errorRowsXenon.append(xenonRow)

    # Test the A matrix
    iodineRowPrecursorSourceSetToZero = matrixPrecursorSourceSetToZero[iodineRow,:]
    iodineRowBaseCase = matrixBaseCase[iodineRow,:]
    # checks the row 
    for columnIndex, val in enumerate(iodineRowPrecursorSourceSetToZero):
        testVal = iodineRowBaseCase[columnIndex]
        if testVal != val:
            errorRowsIodine.append(iodineRow)

    # test xenon b vector
    xenonValPrecursorSourceSetToZero = bVectorPrecursorSourceSetToZero[xenonRow]
    xenonValBaseCase = bVectorBaseCase[xenonRow]
    if xenonValPrecursorSourceSetToZero != xenonValBaseCase:
        errorbVectValsXenon.append(xenonRow)

    # test iodine b vector
    iodineValPrecursorSourceSetToZero = bVectorPrecursorSourceSetToZero[iodineRow]
    iodineValBaseCase = bVectorBaseCase[iodineRow]
    if iodineValPrecursorSourceSetToZero != iodineValBaseCase:
        errorbVectValsIodine.append(iodineRow)

    # test xenon sol vector
    xenonSolValPrecursorSourceSetToZero = solVectorPrecursorSourceSetToZero[xenonRow]
    xenonSolValBaseCase = solVectorBaseCase[xenonRow]
    if xenonSolValPrecursorSourceSetToZero != xenonSolValBaseCase:
        errorSolVectValsXenon.append(xenonRow)

    # test iodine sol vector
    iodineSolValPrecursorSourceSetToZero = solVectorPrecursorSourceSetToZero[iodineRow]
    iodineSolValBaseCase = solVectorBaseCase[iodineRow]
    if iodineSolValPrecursorSourceSetToZero != iodineSolValBaseCase:
        errorSolVectValsIodine.append(iodineRow)
    # Compare xenon soluiton for base case to scipy solution for precursor source set to zero
    compareXeSolutionForBaseCaseToScipySolutionForPrecursorSourceSetToZero.append(
        np.isclose(scipySolutionPrecursorSourceSetToZero[xenonRow], xenonSolValBaseCase))
    # Compare xenon soluiton for base case to scipy solution for base case
    compareXeSolutionForBaseCaseToScipySolutionForBaseCase.append(
            np.isclose(scipySolutionBaseCase[xenonRow], xenonSolValBaseCase))
    # Compare xenon soluiton for precursor source set to zero to scipy solution for base case
    compareXeSolutionForPrecursorSourceSetToZeroToScipySolutionForBaseCase.append(
            np.isclose(scipySolutionBaseCase[xenonRow], xenonSolValPrecursorSourceSetToZero))
    # Compare xenon soluiton for precursor source set to zero to same scipy solution
    compareXeSolutionForPrecursorSourceSetToZeroToSameScipySolution.append(
            np.isclose(scipySolutionPrecursorSourceSetToZero[xenonRow], xenonSolValPrecursorSourceSetToZero))
    # Compare iodine soluiton for base case to scipy solution for precursor source set to zero
    compareISolutionForBaseCaseToScipySolutionForPrecursorSourceSetToZero.append(
            np.isclose(scipySolutionPrecursorSourceSetToZero[iodineRow], iodineSolValBaseCase))
    # Compare iodine soluiton for base case to scipy solution for base case
    compareISolutionForBaseCaseToScipySolutionForBaseCase.append(
            np.isclose(scipySolutionBaseCase[iodineRow], iodineSolValBaseCase))
    # Compare iodine soluiton for precurosr source set to zero to scipy solution for base case
    compareISolutionForPrecursorSourceSetToZeroToScipySolutionForBaseCase.append(
            np.isclose(scipySolutionBaseCase[iodineRow], iodineSolValPrecursorSourceSetToZero))
    # Compare iodine soluiton for precurosr source set to zero to same scipy solution
    compareISolutionForPrecursorSourceSetToZeroToSameScipySolution.append(
            np.isclose(scipySolutionPrecursorSourceSetToZero[iodineRow], iodineSolValPrecursorSourceSetToZero))
    # Compare xenon soluiton for both scipy soluitons
    compareScipySolutions.append(np.isclose(scipySolutionPrecursorSourceSetToZero[xenonRow], scipySolutionBaseCase[xenonRow]))
    # Compare xenon soluiton for both scipy soluiton of only Xe/I
    compareXeIScipySolution.append(np.isclose(scipySolutionPrecursorSourceSetToZero[xenonRow], scipySolutionXeI[k+1]))
    # Compare iodine soluiton for both scipy soluitons
    compareScipySolutions.append(np.isclose(scipySolutionPrecursorSourceSetToZero[iodineRow], scipySolutionBaseCase[iodineRow]))
    # Compare iodine soluiton for both scipy soluiton of only Xe/I
    compareXeIScipySolution.append(np.isclose(scipySolutionPrecursorSourceSetToZero[iodineRow], scipySolutionXeI[k]))
    # Compare xenon solutions from CTF
    compareXenonSolution.append(np.isclose(xenonSolValBaseCase,xenonSolValPrecursorSourceSetToZero))
    # Compare iodine solutions from CTF
    compareIodineSolution.append(np.isclose(iodineSolValBaseCase,iodineSolValPrecursorSourceSetToZero))

    k += 2
"""
