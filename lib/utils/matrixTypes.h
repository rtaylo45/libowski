//*****************************************************************************
// Author: Zack Taylor
//
// Internal type definations used in libowski for matrix stuff.
//*****************************************************************************
#ifndef UTILMATRIXTYPES_H
#define UTILMATRIXTYPES_H
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <complex>

using namespace Eigen;

// Dense matrix stuff
// Complex long double dynamic matrix
typedef Matrix<std::complex<long double>, Dynamic, Dynamic> MatrixCLD;
// Double dynamic matrix
typedef Matrix<double, Dynamic, Dynamic> MatrixD;
// Long double dynamic matrix
typedef Matrix<long double, Dynamic, Dynamic> MatrixLD;

// Useful type definitions for sparce matrix
// Complex long double dynamic matrix
typedef SparseMatrix<std::complex<long double>> SparseMatrixCLD;
// Complex double dynamic matrix
typedef SparseMatrix<std::complex<double>> SparseMatrixCD;
// Double dynamic matrix
typedef SparseMatrix<double> SparseMatrixD;
// Long double dynamic matrix
typedef SparseMatrix<long double> SparseMatrixLD;

#endif
