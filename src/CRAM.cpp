#include "CRAM.h"
#define MTAG1 1

//*************************************************************************
// Initialization of solver
//*************************************************************************
SolverType::SolverType(){
	// Method 0 == CRAM
	// Method 1 == Hyperbolic
	// Method 2 == Parabilic
	int order = 16;
	int method = 0;
	MatrixCLD theta_m16(8,1);
	MatrixCLD alpha_m16(8,1);

	// Defines the complex values for CRAM of order 16
	std::complex<long double> tCRAM1(-10.843917078696988026L, 19.277446167181652284L);
	std::complex<long double> tCRAM2(-5.2649713434426468895L, 16.220221473167927305L);
	std::complex<long double> tCRAM3(5.9481522689511774808L, 3.5874573620183222829L);
	std::complex<long double> tCRAM4(3.5091036084149180974L, 8.4361989858843750826L);
	std::complex<long double> tCRAM5(6.4161776990994341923L, 1.1941223933701386874L);
	std::complex<long double> tCRAM6(1.4193758971856659786L, 10.925363484496722585L);
	std::complex<long double> tCRAM7(4.9931747377179963991L, 5.9968817136039422260L);
	std::complex<long double> tCRAM8(-1.4139284624888862114L, 13.497725698892745389L);


	// Sets the values of theta
	// CRAM
	if (method == 0){
		theta_m16(0,0) = tCRAM1; theta_m16(1,0) = tCRAM2; theta_m16(2,0) = tCRAM3; 
		theta_m16(3,0) = tCRAM4; theta_m16(4,0) = tCRAM5; theta_m16(5,0) = tCRAM6; 
		theta_m16(6,0) = tCRAM7; theta_m16(7,0) = tCRAM8;
		theta = theta_m16;
	}
	// Quad hyperbolic
	else if (method == 1){
		MatrixCLD coeffs = hyperbolicContourCoeffs(order); 
		theta = coeffs.col(0);
	}
	// Quad parabolic
	else if (method == 2){
		MatrixCLD coeffs = parabolicContourCoeffs(order);
		theta = coeffs.col(0);
	}

	// Defines the complex values for CRAM of order 16
	std::complex<long double> aCRAM1(-.0000005090152186522491565L,-.00002422001765285228797L);
	std::complex<long double> aCRAM2(.00021151742182466030907L, .0043892969647380673918L);
	std::complex<long double> aCRAM3(113.39775178483930527L, 101.9472170421585645L);
	std::complex<long double> aCRAM4(15.059585270023467528L, -5.7514052776421819979L);
	std::complex<long double> aCRAM5(-64.500878025539646595L, -224.59440762652096056L);
	std::complex<long double> aCRAM6(-1.4793007113557999718L, 1.7686588323782937906L);
	std::complex<long double> aCRAM7(-62.518392463207918892L, -11.19039109428322848L);
	std::complex<long double> aCRAM8(.041023136835410021273L, -.15743466173455468191L);

	// Sets the values of alpha
	if (method == 0){
		alpha_m16(0,0) = aCRAM1; alpha_m16(1,0) = aCRAM2; alpha_m16(2,0) = aCRAM3; 
		alpha_m16(3,0) = aCRAM4; alpha_m16(4,0) = aCRAM5; alpha_m16(5,0) = aCRAM6; 
		alpha_m16(6,0) = aCRAM7; alpha_m16(7,0) = aCRAM8; 
		alpha_0 = 2.1248537104952237488e-16L;
		alpha = alpha_m16;
	}
	else if (method == 1){
		MatrixCLD coeffs = hyperbolicContourCoeffs(order); 
		alpha = coeffs.col(1);
	}
	else if (method == 2){
		MatrixCLD coeffs = parabolicContourCoeffs(order);
		alpha = coeffs.col(1);
	}

}

//*****************************************************************************
// Calcualtes the quadrature points for the parabolic contour
//
// @param order	Order of the quadrature approx
//*****************************************************************************
MatrixCLD SolverType::parabolicContourCoeffs(int N){
	const std::complex<long double> J(0.0L, 1.0L);
	MatrixCLD coeffs(N/2,2);
	ArrayXLD theta(N/2,1);
	int indexCounter = 0;
	for (int i = 1; i < N; i=i+2){
		theta(indexCounter) = (long double)M_PI*(long double)i/(long double)N;
		indexCounter ++;
	}

	ArrayXCLD phi = (long double)N*(0.1309L - 0.1194L*theta.pow(2) + 0.2500L*theta*J);
	ArrayXCLD phiPrime = (long double)N*(-2.L*0.1194L*theta + 0.2500L*J);
	ArrayXCLD alpha = (J/(long double)N)*phi.exp()*phiPrime;

	// Loops over results to build the return matrix
	for (int i = 0; i < N/2; i++){
		coeffs(i,0) = phi(i);	
		coeffs(i,1) = alpha(i);	
	}
	return coeffs;
}

//*****************************************************************************
// Calcualtes the quadrature points for the hyperbolic contour
//
// @param order	Order of the quadrature approx
//*****************************************************************************
MatrixCLD SolverType::hyperbolicContourCoeffs(int N){
	const std::complex<long double> J(0.0L, 1.0L);
	MatrixCLD coeffs(N/2,2);
	ArrayXLD theta(N/2,1);
	int indexCounter = 0;
	for (int i = 1; i < N; i=i+2){
		theta(indexCounter) = (long double)M_PI*(long double)i/(long double)N;
		indexCounter ++;
	}

	ArrayXCLD phi = (long double)N*2.246L*(1.L - sin(1.1721L - 0.3443L*J*theta));
	ArrayXCLD phiPrime = ((long double)N*cos((theta*3443.L*J)/10000.L - 
		11721.L/10000.L)*3866489.L*J)/5000000.L;
	ArrayXCLD alpha = (J/(long double)N)*phi.exp()*phiPrime;

	// Loops over results to build the return matrix
	for (int i = 0; i < N/2; i++){
		coeffs(i,0) = phi(i);	
		coeffs(i,1) = alpha(i);	
	}
	return coeffs;
}

//*****************************************************************************
// Over rides the solver type 
//
// @param solveType	Type of solve  to use
//							"CRAM"
//							"CRAMScaled"
//*****************************************************************************
void SolverType::setSolveType(std::string solveType){
	assert(solveType == "CRAM" or solveType == "CRAMScaled");

	if (solveType == "CRAM"){
		solverPtr = &SolverType::solveBase;
	}
	else if (solveType == "CRAMScaled"){
		solverPtr = &SolverType::solveScale;
	}
}

//*****************************************************************************
// Matrix expotental solver
//
//	@param A		The coefficient matrix for the system of ODEs
//	@param w0		Initial condition 
//	@param t		Time of the solve
//	
//	return w	Solution vector
//*****************************************************************************
VectorD SolverType::solve(const SparseMatrixD& A, const VectorD& w0, double t){

	return (this->*solverPtr)(A, w0, t);
}
//*****************************************************************************
// Base Matrix expotental solver. No matrix scaling
//	param A		The coefficient matrix for the system of ODEs
//	param w0		Initial condition 
//	param t		Time of the solve
//	
//	return w	Solution vector
//*****************************************************************************
VectorD SolverType::solveBase(const SparseMatrixD& A, const VectorD& w0, double t){

	// The sparse LU solver object
	Eigen::SparseLU<SparseMatrixCLD, COLAMDOrdering<int> > solver;
	
	// MPI stuff
	int myid = mpi.rank;
	int numprocs = mpi.size;
	int eleCount = A.rows();

	// Number of poles
	int s = theta.rows();
	SparseMatrixCLD At(A.rows(),A.cols());
	SparseMatrixCLD tempA(A.rows(),A.cols()); 
	VectorCLD w, tempB, w0cd, myW; 
	VectorD solutionVector;
	w0cd = w0.cast<std::complex<long double>>();
	SparseMatrixCLD ident(A.rows(),A.cols());

	myW = 0.*w0cd, w = 0*w0cd;
	ident.setIdentity();
	At = A.cast<std::complex<long double>>()*t;
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
	solutionVector = w.real().cast<double>();

	return solutionVector;
}

//*****************************************************************************
// Scaled Matrix expotental solver.
//	param A		The coefficient matrix for the system of ODEs
//	param w0		Initial condition 
//	param t		Time of the solve
//	
//	return w	Solution vector
//*****************************************************************************
VectorD SolverType::solveScale(const SparseMatrixD& A, const VectorD& w0, double t){

	// MPI stuff
	int myid = mpi.rank;
	int numprocs = mpi.size;
	int eleCount = A.rows()*A.rows();
	int beta = 10;
	double matrixReduction = pow(2.,beta);
	int matrixPower = pow(2,beta);

	// Number of poles
	int s = theta.rows();
	SparseMatrixCLD At(A.rows(),A.cols());
	SparseMatrixCLD tempA(A.rows(),A.cols()); 
	SparseMatrixCLD myExpA;
	SparseMatrixCLD expA;
	SparseMatrixCLD ident(A.rows(),A.cols());
	SparseMatrixLD expAReal;
	SparseMatrixLD expASquaredReal;
	VectorD solutionVector;

	myExpA = 0*A.cast<std::complex<long double>>();
	ident.setIdentity();
	expA = 0*A.cast<std::complex<long double>>();
	At = A.cast<std::complex<long double>>()*t/matrixReduction;

	// Loops over the imaginary poles to compute the matrix exponential
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
	solutionVector = expASquaredReal.cast<double>()*w0.cast<double>();

	return solutionVector;
}
