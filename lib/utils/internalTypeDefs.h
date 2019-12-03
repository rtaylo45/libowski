//*****************************************************************************
// Author: Zack Taylor
//
// Internal type definations used in libowski that eigen doesn't define
//*****************************************************************************
#ifndef UTILINTERNALTYPEDEFS_H
#define UTILINTERNALTYPEDEFS_H
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <complex>

using namespace Eigen;

// Useful type definitions for dense matrices and vectors
// Complex long double dynamic matrix
typedef Matrix<std::complex<long double>, Dynamic, Dynamic> MatrixXcld;
// Complex long double dynamic vector
typedef Matrix<std::complex<long double>, Dynamic, Dynamic> VectorXcld;

// Useful type definitions for sparce matrices and vectors
// Complex long double dynamic matrix
typedef SparceMatrix<std::complex<long double>> SpaceMatrixXcld;
// Complex double dynamic matrix
typedef SparceMatrix<std::complex<double>> SparceMatrixXcd;
// Double dynamic matrix
typedef SparceMatrix<double> SparceMatrixXd;
// Long double dynamic matrix
typedef SparceMatrix<long double> SparceMatrixXld;
// double dynamic vector
typedef SparceVector<double> SparceVectorXd;
// long double dynamic vector
typedef SparceVector<long double> SparceVectorXld;


#endif
