//*****************************************************************************
// Author: Zack Taylor
//
// Linear algebra utilities that Eigen doesn't provide
//*****************************************************************************
#include "linearAlgebra.h"

//*****************************************************************************
// Implements the Moore-Penrose pseudo inverse for a eigen sparse matrix
//
// @param A		matrix
//*****************************************************************************
template <typename derived>
SparseMatrix<derived> MoorePenroseInv(SparseMatrix<derived> A){
	SparseMatrix<derived> Ainv;
	SparseMatrix<derived> AConjugateTranspose;
	SparseMatrix<derived> Atemp;
	SparseMatrix<derived> I(A.rows(),A.cols());
	SparseQR<SparseMatrix<derived>, COLAMDOrdering<int>> solver;

	I.setIdentity();	
	AConjugateTranspose = A.adjoint();
	Atemp = AConjugateTranspose*A;
	solver.compute(Atemp);
	Ainv = solver.solve(I)*AConjugateTranspose;
	return Ainv;
}

//*****************************************************************************
// Matrix Squaring for sparse matrix
//
// @param A			Matrix to be squared
// @param alpha	Matrix power
//*****************************************************************************
template <typename derived>
SparseMatrix<derived> MatrixSquare(SparseMatrix<derived> A, int alpha){
	SparseMatrix<derived> Areturn;
	SparseMatrix<derived> ASquared;
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

// Data types that can use the template functions
template SparseMatrixLD MoorePenroseInv(SparseMatrixLD A);
template SparseMatrixD MoorePenroseInv(SparseMatrixD A);
template SparseMatrixCLD MoorePenroseInv(SparseMatrixCLD A);
template SparseMatrixCD MoorePenroseInv(SparseMatrixCD A);

template SparseMatrixLD MatrixSquare(SparseMatrixLD A, int alpha);
template SparseMatrixD MatrixSquare(SparseMatrixD A, int alpha);
