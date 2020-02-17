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


//*****************************************************************************
// Pseudo inverse linearly independent columns for sparse matrix
//*****************************************************************************
template <typename derived>
SparseMatrix<derived> MoorePenroseInv(const SparseMatrix<derived>&);

//*****************************************************************************
// Binomial Coefficient
//*****************************************************************************
int binomialCoeff(int, int);

//*****************************************************************************
// Factorial
//*****************************************************************************
int factorial(int);

//*****************************************************************************
// l1 norm sparse matrix
//*****************************************************************************
double l1norm(const SparseMatrixD&);

//*****************************************************************************
// Arnoldi algorithm
//*****************************************************************************
void arnoldi(const SparseMatrixD& A, const VectorD& b, const int n, 
	MatrixD& V, SparseMatrixD& H);

#endif
