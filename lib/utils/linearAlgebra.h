//*****************************************************************************
// Author: Zack Taylor
//
// Linear algebra utilities that Eigen doesn't provide
//*****************************************************************************
#ifndef UTILLINEARALGEBRA_H
#define UTILLINEARALGEBRA_H
#include "matrixTypes.h"
#include "vectorTypes.h"

using namespace Eigen;

// Pseudo inverse linearly independent columns for complex long double sparse matrix
SparseMatrixCLD MoorePenroseInv(SparseMatrixCLD);
// Pseudo inverse linearly independent columns for complex double sparse matrix
SparseMatrixCD MoorePenroseInv(SparseMatrixCD);
// Pseudo inverse lineary independent columns for long double sparse matrix
SparseMatrixLD MoorePenroseInv(SparseMatrixLD);
// Pseudo inverse lineary independent columns for double sparse matrix
SparseMatrixD MoorePenroseInv(SparseMatrixD);

// Compute matrix squaring for long double sparse matrix
SparseMatrixLD MatrixSquare(SparseMatrixLD, int);
// Compute matrix squaring for double sparse matrix
SparseMatrixD MatrixSquare(SparseMatrixD, int);
#endif
