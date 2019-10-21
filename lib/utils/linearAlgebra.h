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

// Pseudo inverse for complex double sparse matrix
SparseMatrix<std::complex<double>> MoorePenroseInv(
	SparseMatrix<std::complex<double>>);

// Pseudo inverse for double sparse matrix
SparseMatrix<double> MoorePenroseInv(SparseMatrix<double>);
#endif
