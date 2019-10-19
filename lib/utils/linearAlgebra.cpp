//*****************************************************************************
// Author: Zack Taylor
//
// Linear algebra utilities that Eigen doesn't provide
//*****************************************************************************
#include "linearAlgebra.h"

//*****************************************************************************
// Implements the Moore-Penrose pseudo inverse for a eigen sparse 
// complex matrix
//
// @param A		Complex matrix
//*****************************************************************************
SparseMatrix<std::complex<double>> MoorePenroseInv(
	SparseMatrix<std::complex<double>> A){
	SparseMatrix<std::complex<double>> Ainv;
	SparseMatrix<std::complex<double>> AConjugateTranspose;
	SparseMatrix<std::complex<double>> Atemp;
	SparseMatrix<std::complex<double>> I(A.rows(),A.cols());
	SparseLU<SparseMatrix<std::complex<double>>, COLAMDOrdering<int>> solver;

	I.setIdentity();	
	AConjugateTranspose = A.adjoint();
	Atemp = A*AConjugateTranspose;
	solver.factorize(Atemp);
	Ainv = AConjugateTranspose*solver.solve(I);
	return Ainv;
}
