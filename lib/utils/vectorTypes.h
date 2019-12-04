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
// Complex long double dynamic vector
typedef Matrix<std::complex<long double>, Dynamic, 1> VectorCLD;
// Double dynamic vector
typedef Matrix<double, Dynamic, 1> VectorD;
// Long double dynamic vector
typedef Matrix<long double, Dynamic, 1> VectorLD;

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
