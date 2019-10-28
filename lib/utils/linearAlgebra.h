//*****************************************************************************
// Author: Zack Taylor
//
// Linear algebra utilities that Eigen doesn't provide
//*****************************************************************************
#ifndef UTILLINEARALGEBRA_H
#define UTILLINEARALGEBRA_H
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <complex>
#include <iostream>

using namespace Eigen;

// Pseudo inverse linearly independent columns for complex double sparse matrix
SparseMatrix<std::complex<double>> MoorePenroseInv(SparseMatrix<std::complex<double>>);
// Pseudo inverse lineary independent columns for double sparse matrix
SparseMatrix<double> MoorePenroseInv(SparseMatrix<double>);

// Compute matrix squaring for complex double sparse matrix
//SparseMatrix<std::complex<double>> MatrixSquare(SparseMatrix<std::complex<double>>, int);
// Compute matrix squaring for double sparse matrix
SparseMatrix<double> MatrixSquare(SparseMatrix<double>, int);
#endif
