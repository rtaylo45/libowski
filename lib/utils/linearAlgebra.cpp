//*****************************************************************************
// Author: Zack Taylor
//
// Linear algebra utilities that Eigen doesn't provide
//*****************************************************************************
#include "linearAlgebra.h"

//*****************************************************************************
// Implements the Moore-Penrose pseudo inverse for a eigen sparse 
// complex long double matrix
//
// @param A		Complex matrix
//*****************************************************************************
SparseMatrixCLD MoorePenroseInv(SparseMatrixCLD A){
	SparseMatrixCLD Ainv;
	SparseMatrixCLD AConjugateTranspose;
	SparseMatrixCLD Atemp;
	SparseMatrixCLD I(A.rows(),A.cols());
	SparseQR<SparseMatrixCLD, COLAMDOrdering<int>> solver;

	I.setIdentity();	
	AConjugateTranspose = A.adjoint();
	Atemp = AConjugateTranspose*A;
	solver.compute(Atemp);
	Ainv = solver.solve(I)*AConjugateTranspose;
	return Ainv;
}

//*****************************************************************************
// Implements the Moore-Penrose pseudo inverse for a eigen sparse 
// complex double matrix
//
// @param A		Complex matrix
//*****************************************************************************
SparseMatrixCD MoorePenroseInv(SparseMatrixCD A){
	SparseMatrixCD Ainv;
	SparseMatrixCD AConjugateTranspose;
	SparseMatrixCD Atemp;
	SparseMatrixCD I(A.rows(),A.cols());
	SparseQR<SparseMatrixCD, COLAMDOrdering<int>> solver;

	I.setIdentity();	
	AConjugateTranspose = A.adjoint();
	Atemp = AConjugateTranspose*A;
	solver.compute(Atemp);
	Ainv = solver.solve(I)*AConjugateTranspose;
	return Ainv;
}

//*****************************************************************************
// Implements the Moore-Penrose pseudo inverse for a eigen sparse 
// long double matrix
//
// @param A		double sparse matrix
//*****************************************************************************
SparseMatrixLD MoorePenroseInv(SparseMatrixLD A){
	SparseMatrixLD Ainv;
	SparseMatrixLD AConjugateTranspose;
	SparseMatrixLD Atemp;
	SparseMatrixLD I(A.rows(),A.cols());
	SparseQR<SparseMatrixLD, COLAMDOrdering<int>> solver;

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
SparseMatrixD MoorePenroseInv(SparseMatrixD A){
	SparseMatrixD Ainv;
	SparseMatrixD AConjugateTranspose;
	SparseMatrixD Atemp;
	SparseMatrixD I(A.rows(),A.cols());
	SparseQR<SparseMatrixD, COLAMDOrdering<int>> solver;

	I.setIdentity();	
	AConjugateTranspose = A.adjoint();
	Atemp = AConjugateTranspose*A;
	solver.compute(Atemp);
	Ainv = solver.solve(I)*AConjugateTranspose;
	return Ainv;
}

//*****************************************************************************
// Matrix Squaring for sparse long double matrix
//
// @param A			Matrix to be squared
// @param alpha	Matrix power
//*****************************************************************************
SparseMatrixLD MatrixSquare(SparseMatrixLD A, int alpha){
	SparseMatrixLD Areturn;
	SparseMatrixLD ASquared;
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

//*****************************************************************************
// Matrix Squaring for sparse double matrix
//
// @param A			Matrix to be squared
// @param alpha	Matrix power
//*****************************************************************************
SparseMatrixD MatrixSquare(SparseMatrixD A, int alpha){
	SparseMatrixD Areturn;
	SparseMatrixD ASquared;
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
