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
	MatrixXcd w, tempB, w0cd, myW; 
	w0cd = w0.cast<std::complex<double>>();
	SparseMatrix<double> ident = buildSparseIdentity(A.rows());

	// MPI stuff
	int myid, numprocs, islave;
	MPI_Status status;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	std::cout << numprocs << myid << std::endl;

	w = 0.*w0cd;
	At = A.cast<std::complex<double>>()*t;

	//for (;;){
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
			MPI_Send(&myW, 1, MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
		}
		else {
			w = myW;
			for (islave = 1; islave < numprocs; islave++) {
				MPI_Recv(&myW, 1, MPI_DOUBLE_COMPLEX, islave, 1, MPI_COMM_WORLD, &status);
				w += myW;
			}
		}
	//}
	MPI_Finalize();
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
