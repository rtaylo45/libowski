#include "coreSolver.h"
#include <iostream>

using namespace Eigen;

//*****************************************************************************
// Matrix expotental solver
//	param A		The coefficient matrix for the system of ODEs
//	param w0	Initial condition 
//	param t		Time of the solve
//	
//	return w	Solution vector
//*****************************************************************************
MatrixXd SolverType::solve(SparseMatrix<double> A, VectorXcd w0, double t){
	// The sparse LU solver object
	SparseLU<SparseMatrix<std::complex<double>>, COLAMDOrdering<int> > solver;

	// Number of poles
	int s = 8;
	SparseMatrix<std::complex<double>> At(A.rows(),A.cols()); 
	SparseMatrix<std::complex<double>> tempA(A.rows(),A.cols()); 
	VectorXcd w, tempB, w0cd; 
	w0cd = w0.cast<std::complex<double>>();
	At = A.cast<std::complex<double>>()*t;
	w = 0.*w0cd;
	SparseMatrix<double> ident = buildSparseIdentity(A.rows());

	for (int k = 0; k < s; k++){
		tempA = At - theta(k)*ident;
		tempB = alpha(k)*w0cd;
		// analyze the sparsisty pattern
		solver.analyzePattern(tempA);
		// Compute the numerical factorization
		solver.factorize(tempA);

		w = w + solver.solve(tempB);
	}
	w = 2.*w.real();
	w = w + alpha_0*w0cd;

	return w.real();
}

//*****************************************************************************
// Builds a sparse identity matrix
//	param n		Square matrix size
//
//	return nxn	identity matrix
//*****************************************************************************
SparseMatrix<double> SolverType::buildSparseIdentity(int n){

	SparseMatrix<double> ident(n,n);
	for (int i = 0; i < n; i++){
		ident.insert(i,i) = 1.0;
	}
	return ident;
}
