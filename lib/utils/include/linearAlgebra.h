//*****************************************************************************
// Author: Zack Taylor
//
// Linear algebra utilities that Eigen doesn't provide
//*****************************************************************************
#ifndef UTILLINEARALGEBRA_H
#define UTILLINEARALGEBRA_H
#include "matrixTypes.h"
#include "vectorTypes.h"
#include "exception.h"

using namespace Eigen;

//*****************************************************************************
// Pseudo inverse linearly independent columns for sparse matrix
//*****************************************************************************
template <typename derived>
SparseMatrix<derived> MoorePenroseInv(const SparseMatrix<derived>&);


//*****************************************************************************
// l1 norm sparse matrix
//*****************************************************************************
double l1norm(const SparseMatrixD&);

//*****************************************************************************
// Produces the l1norm of A*B. Need to eventually change this to use the
// estiment of the norm which is use in the paper.
//*****************************************************************************
double normAm(const SparseMatrixD&, const SparseMatrixD&);

//*****************************************************************************
// Produces the 1lnorm of A^m
//*****************************************************************************
double normAm(const SparseMatrixD&, const int);

//*****************************************************************************
// Arnoldi algorithm
//*****************************************************************************
void arnoldi(const SparseMatrixD& A, const VectorD& b, const int n,
  MatrixD& V, SparseMatrixD& H);

#endif
