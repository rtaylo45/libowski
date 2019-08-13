#include "coreSolver.h"
#include <iostream>
#define MTAG1 1

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

	// MPI stuff
	int myid, numprocs, ierr;
	//MPI_Init(NULL, NULL);
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	int eleCount = A.rows();

	// Number of poles
	int s = 8;
	SparseMatrix<std::complex<double>> At(A.rows(),A.cols()); 
	SparseMatrix<std::complex<double>> tempA(A.rows(),A.cols()); 
	VectorXcd w, tempB, w0cd, myW; 
	w0cd = w0.cast<std::complex<double>>();
	SparseMatrix<double> ident = buildSparseIdentity(A.rows());


	myW = 0.*w0cd, w = 0*w0cd;
	At = A.cast<std::complex<double>>()*t;

	for (int k = myid; k < s; k += numprocs){
		tempA = At - theta(k)*ident;
		tempB = alpha(k)*w0cd;
		// analyze the sparsisty pattern
		solver.analyzePattern(tempA);
		// Compute the numerical factorization
		solver.factorize(tempA);

		myW = myW + solver.solve(tempB);
	}
	if (myid != 0){
		ierr = MPI_Send(myW.data(), eleCount, MPI::DOUBLE_COMPLEX, 0, MTAG1, MPI_COMM_WORLD);
	}
	else {
		w = myW;
		for (int islave = 1; islave < numprocs; islave++) {
			ierr = MPI_Recv(myW.data(), eleCount, MPI::DOUBLE_COMPLEX, islave, MTAG1, MPI_COMM_WORLD, &status);
			w = w + myW;
		}
	}
	//MPI_Finalize();
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
