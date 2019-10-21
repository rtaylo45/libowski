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
	SparseQR<SparseMatrix<std::complex<double>>, COLAMDOrdering<int>> solver;

	I.setIdentity();	
	AConjugateTranspose = A.adjoint();
	Atemp = A*AConjugateTranspose;
	solver.compute(Atemp);
	Ainv = AConjugateTranspose*solver.solve(I);
	return Ainv;
}

//*****************************************************************************
// Implements the Moore-Penrose pseudo inverse for a eigen sparse 
// double matrix
//
// @param A		double sparse matrix
//*****************************************************************************
SparseMatrix<double> MoorePenroseInv(SparseMatrix<double> A){
	SparseMatrix<double> Ainv;
	SparseMatrix<double> AConjugateTranspose;
	SparseMatrix<double> Atemp;
	SparseMatrix<double> I(A.rows(),A.cols());
	SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;

	I.setIdentity();	
	AConjugateTranspose = A.adjoint();
	Atemp = A*AConjugateTranspose;
	solver.analyzePattern(Atemp);
	solver.factorize(Atemp);
	//solver.compute(Atemp);
	Ainv = AConjugateTranspose*solver.solve(I);
	return Ainv;
}
