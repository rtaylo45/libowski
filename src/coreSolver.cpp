#include "coreSolver.h"
#include <iostream>

using namespace Eigen;

// Solver
MatrixXd SolverType::solve(SparseMatrix<double> A, VectorXcd w0, double t){
	SparseLU<SparseMatrix<std::complex<double>>, COLAMDOrdering<int> >   solver;

	int s = 8;
	SparseMatrix<std::complex<double>> At(A.rows(),A.cols()); 
	SparseMatrix<std::complex<double>> tempA(A.rows(),A.cols()); 
	VectorXcd w, tempB, w0cd; 
	w0cd = w0.cast<std::complex<double>>();
	At = A.cast<std::complex<double>>()*t;
	w = 0.*w0cd;
	SparseMatrix<double> ident = buildSparseIdentity(A.rows());
	// Compute the ordering permutation vector from the structural pattern of A

	for (int k = 0; k < s; k++){
		tempA = At - theta(k)*ident;
		tempB = alpha(k)*w0cd;
		// Compute the numerical factorization
		solver.analyzePattern(tempA);
		solver.factorize(tempA);

		w = w + solver.solve(tempB);
	}
	w = 2.*w.real();
	w = w + alpha_0*w0cd;

	return w.real();
}

// Builds a sparse identity matrix
SparseMatrix<double> SolverType::buildSparseIdentity(int n){

	SparseMatrix<double> ident(n,n);
	for (int i = 0; i < n; i++){
		ident.insert(i,i) = 1.0;
	}
	return ident;
}
