#include "CRAM.h"
#define MTAG1 1

using namespace Eigen;
//*************************************************************************
// Initialization of solver
//*************************************************************************
SolverType::SolverType(){
	// Defines the complex values
	std::complex<double> t1(-10.843917078696988026, 19.277446167181652284);
	std::complex<double> t2(-5.2649713434426468895, 16.220221473167927305);
	std::complex<double> t3(5.9481522689511774808, 3.5874573620183222829);
	std::complex<double> t4(3.5091036084149180974, 8.4361989858843750826);
	std::complex<double> t5(6.4161776990994341923, 1.1941223933701386874);
	std::complex<double> t6(1.4193758971856659786, 10.925363484496722585);
	std::complex<double> t7(4.9931747377179963991, 5.9968817136039422260);
	std::complex<double> t8(-1.4139284624888862114, 13.497725698892745389);

	// Sets the values of theta
	theta(0,0) = t1; theta(1,0) = t2; theta(2,0) = t3; theta(3,0) = t4;
	theta(4,0) = t5; theta(5,0) = t6; theta(6,0) = t7; theta(7,0) = t8;

	// Defines the complex values
	std::complex<double> a1(-.0000005090152186522491565,-.00002422001765285228797);
	std::complex<double> a2(.00021151742182466030907, .0043892969647380673918);
	std::complex<double> a3(113.39775178483930527, 101.9472170421585645);
	std::complex<double> a4(15.059585270023467528, -5.7514052776421819979);
	std::complex<double> a5(-64.500878025539646595, -224.59440762652096056);
	std::complex<double> a6(-1.4793007113557999718, 1.7686588323782937906);
	std::complex<double> a7(-62.518392463207918892, -11.19039109428322848);
	std::complex<double> a8(.041023136835410021273, -.15743466173455468191);

	// Sets the values of alpha
	alpha(0,0) = a1; alpha(1,0) = a2; alpha(2,0) = a3; alpha(3,0) = a4; 
	alpha(4,0) = a5; alpha(5,0) = a6; alpha(6,0) = a7; alpha(7,0) = a8; 
}

//*****************************************************************************
// Matrix expotental solver
//	param A		The coefficient matrix for the system of ODEs
//	param w0		Initial condition 
//	param t		Time of the solve
//	
//	return w	Solution vector
//*****************************************************************************
MatrixXd SolverType::solve(SparseMatrix<double> A, VectorXd w0, double t){

	return solveBase(A, w0, t);
}
//*****************************************************************************
// Base Matrix expotental solver. No matrix scaling
//	param A		The coefficient matrix for the system of ODEs
//	param w0		Initial condition 
//	param t		Time of the solve
//	
//	return w	Solution vector
//*****************************************************************************
MatrixXd SolverType::solveBase(SparseMatrix<double> A, VectorXd w0, double t){

	// The sparse LU solver object
	SparseLU<SparseMatrix<std::complex<double>>, COLAMDOrdering<int> > solver;

	// MPI stuff
	int myid = mpi.rank;
	int numprocs = mpi.size;
	int eleCount = A.rows();

	// Number of poles
	int s = 8;
	SparseMatrix<std::complex<double>> At(A.rows(),A.cols());
	SparseMatrix<std::complex<double>> tempA(A.rows(),A.cols()); 
	VectorXcd w, tempB, w0cd, myW; 
	w0cd = w0.cast<std::complex<double>>();
	SparseMatrix<std::complex<double>> ident = buildSparseIdentity(A.rows());

	myW = 0.*w0cd, w = 0*w0cd;
	At = A.cast<std::complex<double>>()*t;
	// analyze the sparsisty pattern
	solver.analyzePattern(tempA);

	// Loops over the imaginary poles. This is a linear solve over 8 lineary 
	// independent systems. The sum of all the independent solutions is w.
	for (int k = myid; k < s; k += numprocs){
		tempA = At - theta(k)*ident;
		tempB = alpha(k)*w0cd;
		// Compute the numerical factorization
		solver.factorize(tempA);

		myW = myW + solver.solve(tempB);
	}
	if (myid != 0){
		// Sends solution data to the master node 
		mpi.send(myW, eleCount, 0, MTAG1);
	}
	else {
		w = myW;
		// Receives data from the slave nodes
		for (int islave = 1; islave < numprocs; islave++) {
			myW = mpi.recv(myW, eleCount, islave, MTAG1);
			w = w + myW;
		}
	}
	w = 2.*w.real();
	w = w + alpha_0*w0cd;

	return w.real();
}

//*****************************************************************************
// Scaled Matrix expotental solver.
//	param A		The coefficient matrix for the system of ODEs
//	param w0		Initial condition 
//	param t		Time of the solve
//	
//	return w	Solution vector
//*****************************************************************************
MatrixXd SolverType::solveScale(SparseMatrix<double> A, VectorXd w0, double t){

	// MPI stuff
	int myid = mpi.rank;
	int numprocs = mpi.size;
	int eleCount = A.rows()*A.rows();
	int beta = 10;
	double matrixReduction = pow(2.,beta);
	int matrixPower = pow(2,beta);

	// Number of poles
	int s = 8;
	SparseMatrix<std::complex<double>> At(A.rows(),A.cols());
	SparseMatrix<std::complex<double>> tempA(A.rows(),A.cols()); 
	SparseMatrix<std::complex<double>> myExpA;
	SparseMatrix<std::complex<double>> expA;
	SparseMatrix<std::complex<double>> ident = buildSparseIdentity(A.rows());
	SparseMatrix<double> expAReal;
	SparseMatrix<double> expASquaredReal;

	myExpA = 0*A.cast<std::complex<double>>();
	expA = 0*A.cast<std::complex<double>>();
	At = A.cast<std::complex<double>>()*t/matrixReduction;
	//std::cout << At << " " << matrixReduction << std::endl;
	//std::cout << " " << std::endl;

	// Loops over the imaginary poles. This is a linear solve over 8 lineary 
	// independent systems. The sum of all the independent solutions is w.
	for (int k = myid; k < s; k += numprocs){
		tempA = At - theta(k)*ident;

		// Multipy A-1*alpha(K) because the matrix is 1-1 and only left invertible 
		// because it has linearly independent columns but not rows. 
		myExpA = myExpA + MoorePenroseInv(tempA)*alpha(k);
	}
	if (myid != 0){
		// Sends solution data to the master node 
		mpi.send(myExpA, eleCount, 0, MTAG1);
	}
	else {
		expA = myExpA;
		// Receives data from the slave nodes
		for (int islave = 1; islave < numprocs; islave++) {
			myExpA = mpi.recv(myExpA, eleCount, islave, MTAG1);
			expA = expA + myExpA;
		}
	}
	expAReal = 2.*expA.real();
	expAReal = expAReal + alpha_0*ident.real();
	expASquaredReal = MatrixSquare(expAReal, matrixPower);

	return expASquaredReal*w0;
}
//*****************************************************************************
// Builds a sparse identity matrix
//	param n		Square matrix size
//
//	return nxn	identity matrix
//*****************************************************************************
SparseMatrix<std::complex<double>> SolverType::buildSparseIdentity(int n){

	SparseMatrix<std::complex<double>> ident(n,n);
	for (int i = 0; i < n; i++){
		ident.insert(i,i) = 1.0;
	}
	return ident;
}
