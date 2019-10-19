//*****************************************************************************
// Author: Zack Taylor
//
// Linear algebra utilities that Eigen doesn't provide
//*****************************************************************************
#ifndef UTILLINEARALGEBRA_H
#define UTILLINEARALGEBRA_H
#include <Eigen/Sparse>

using namespace Eigen;

SparseMatrix<std::complex<double>> MoorePenroseInv(
	SparseMatrix<std::complex<double>>);
#endif
