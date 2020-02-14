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
SparseMatrix<derived> MoorePenroseInv(const SparseMatrix<derived>& A){
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
SparseMatrix<derived> MatrixSquare(const SparseMatrix<derived>& A, int alpha){
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

//*****************************************************************************
// Binomial coefficient. N choose k
// Taken from:
// https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
//
// @param N		Number of things
// @param k		Number to be taken
//*****************************************************************************
int binomialCoeff(int n, int k){  
    int res = 1;  
  
   // Since C(n, k) = C(n, n-k)  
   if ( k > n - k ){
      k = n - k;  
	} 
   // Calculate value of  
   // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]  
   for (int i = 0; i < k; ++i){  
      res *= (n - i);  
      res /= (i + 1);  
   }  
  
   return res;  
} 
//*****************************************************************************
// Factorial
//
// @param n		Factorial number
//*****************************************************************************
int factorial(int n)
{
   if(n > 1){
		return n * factorial(n - 1);
	}
   else{
		return 1;
	}		
}

//*****************************************************************************
// l1 norm of a sparse matrix
//
// @param A		Sparse matrix
//*****************************************************************************
double l1norm(const SparseMatrixD& A){
	double norm;
	norm = (A.cwiseAbs()*VectorD::Ones(A.cols())).maxCoeff();
	return norm;

}

// Data types that can use the template functions
template SparseMatrixLD MoorePenroseInv(const SparseMatrixLD& A);
template SparseMatrixD MoorePenroseInv(const SparseMatrixD& A);
template SparseMatrixCLD MoorePenroseInv(const SparseMatrixCLD& A);
template SparseMatrixCD MoorePenroseInv(const SparseMatrixCD& A);

template SparseMatrixLD MatrixSquare(const SparseMatrixLD& A, int alpha);
template SparseMatrixD MatrixSquare(const SparseMatrixD& A, int alpha);

