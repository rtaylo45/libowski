//*****************************************************************************
// Author: Zack Taylor
//
// Internal type definations used in libowski for vector stuff.
//*****************************************************************************
#ifndef UTILVECTORTYPES_H
#define UTILVECTORTYPES_H
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <complex>

using namespace Eigen;

// Dense matrix stuff
// Complex long double dynamic vector
typedef Matrix<std::complex<long double>, Dynamic, 1> VectorCLD;
// Complex double dynamic vector
typedef Matrix<std::complex<double>, Dynamic, 1> VectorCD;
// Double dynamic vector
typedef Matrix<double, Dynamic, 1> VectorD;
// Long double dynamic vector
typedef Matrix<long double, Dynamic, 1> VectorLD;
// Integer vector
typedef Matrix<int, Dynamic, 1> VectorI;
// Long integer vector
typedef Matrix<long int, Dynamic, 1> VectorLI;

// Useful type definitions for sparce vector
// Complex long double dynamic vector
typedef SparseVector<std::complex<long double>> SpaseVectorCLD;
// Complex double dynamic vector
typedef SparseVector<std::complex<double>> SparseVectorCD;
// Double dynamic vector
typedef SparseVector<double> SparseVectorD;
// Long double dynamic vector
typedef SparseVector<long double> SparseVectorLD;

#endif
