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

// Dense Array stuff. Arrays in eigen allow for element wise operations.
// They are useful for holding values that i can manipulate.
// Long double dynamic array
typedef Array<long double, Dynamic, Dynamic> ArrayXLD;
// Complex long double dynamic array
typedef Array<std::complex<long double>, Dynamic, Dynamic> ArrayXCLD;
// Double dynamic array
typedef Array<double, Dynamic, Dynamic> ArrayD;

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
