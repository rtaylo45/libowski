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
	Atemp = AConjugateTranspose*A;
	solver.compute(Atemp);
	Ainv = solver.solve(I)*AConjugateTranspose;
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
	Atemp = AConjugateTranspose*A;
	solver.compute(Atemp);
	Ainv = solver.solve(I)*AConjugateTranspose;
	return Ainv;
}

//*****************************************************************************
// Matrix Squaring for sparse double matrix
//
// @param A			Matrix to be squared
// @param alpha	Matrix power
//*****************************************************************************
SparseMatrix<double> MatrixSquare(SparseMatrix<double> A, int alpha){
	SparseMatrix<double> Areturn;
	SparseMatrix<double> ASquared;
	Areturn = A;
	
	if (alpha != 1){
		ASquared = A*A;
		Areturn = ASquared;
		
		// Loops over the number of times to square the matrix	
		for (int i = 0; i < alpha/2-1; i++){
			Areturn = Areturn*ASquared;	
		}
	}
	return Areturn;
}
