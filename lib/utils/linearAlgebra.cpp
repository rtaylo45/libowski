//*****************************************************************************
// Author: Zack Taylor
//
// Linear algebra utilities that Eigen doesn't provide
//*****************************************************************************
#include "linearAlgebra.h"
#include <iostream>

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

//*****************************************************************************
// Arnoldi algorithm
//
// Used for generating the Krylov subspace. This process produces and 
// orthnormal basis V(m x n) and an upper Hessenberg matrix (n x n). The 
// algorithm produces coefficinet for h(n+1,n) and vector v(n+1) but these
// are not added to the returned matrices. In the future these can be added
// to produced to more accurate Krylov subspace approximation.
//
// The orthnormal basis is a dense matrix and will be retured as one H is 
// lower triangular so it is retured as sparse. H is also the matrix that 
// will be evaluated with the matrix exponential so it needs to be sparse to
// be evaluated by the existing matrix exp methods. 
//
// @param A		Sparse matrix 
// @param b		b vector used to build the subspace
// @param n		Dimention of the subspace
// @param V		Returned orthonormal basis
// @param H		Returned upper Hessenburg matrix, the projection of A into the 
//					Krylov subspace
//*****************************************************************************
void arnoldi(const SparseMatrixD& A, const VectorD& b, const int n, 
	MatrixD& V, SparseMatrixD& H){
	VectorD q, v;
	double h, normV;
	const int m = A.cols();
	const double eps = 1.e-12;

	// Set temp matrices to build the subspace
	MatrixD Vtemp = MatrixD::Zero(m, n+1);
	SparseMatrixD Htemp(n+1, n);
	// reserve space for the H matrix. This includes the main diagonal 
	// the sub main diagonal, and all elements above the main diagonal.
	Htemp.reserve(n*(n+1) - (n-1)*n/2 - 1);

	q = b/b.norm();
	Vtemp.col(0) = q;

	// Loop over columns
	for (int i = 0; i < n; i++){
		// new candidant vector
		v = A*q;
		// Subtract previous vectors
		for (int j = 0; j <= i; j++){
			h = Vtemp.col(j).conjugate().dot(v);
			Htemp.insert(j, i) = h;
			v = v - h*Vtemp.col(j);

		}
		normV = v.norm();
		Htemp.insert(i+1,i) = normV;
		// If v is too small its the zero vector and should not be added
		if (normV > eps){
			q = v/normV;
			Vtemp.col(i+1) = q;
		}
		else{
			V = Vtemp.topLeftCorner(n,n);
			H = Htemp.topLeftCorner(n,n);
			return;
		}
	}
	V = Vtemp.topLeftCorner(n,n);
	H = Htemp.topLeftCorner(n,n);
	return;
}	

// Data types that can use the template functions
template SparseMatrixLD MoorePenroseInv(const SparseMatrixLD& A);
template SparseMatrixD MoorePenroseInv(const SparseMatrixD& A);
template SparseMatrixCLD MoorePenroseInv(const SparseMatrixCLD& A);
template SparseMatrixCD MoorePenroseInv(const SparseMatrixCD& A);

