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

// Pseudo inverse linearly independent columns for sparse matrix
template <typename derived>
SparseMatrix<derived> MoorePenroseInv(SparseMatrix<derived>);

// Compute matrix squaring for sparse matrix
template <typename derived>
SparseMatrix<derived> MatrixSquare(SparseMatrix<derived>, int);
#endif
