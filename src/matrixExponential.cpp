#include "matrixExponential.h"
#define MTAG1 1
//*****************************************************************************
// Methods for MatrixExponential factory class
//*****************************************************************************
matrixExponential *matrixExponentialFactory::getExpSolver(std::string type){
	matrixExponential *solver = nullptr;

	if (type == "CRAM"){
		solver = new CRAM;
		return solver;
	}
	else if (type == "parabolic"){
		solver = new parabolic;
		return solver;
	}
	else if (type == "hyperbolic"){
		solver = new hyperbolic;
		return solver;
	}
	else {
		std::cout << "You fucked up and tried to pick a matrix \n"
			"exponential solver that Zack has not put in. \n"
			"Please constact Zack at 1 800 eat shit" << std::endl;
		exit(1);
	}
}

//*****************************************************************************
// Methods for Cauchy Class
//*****************************************************************************

//*****************************************************************************
// Calculates exp(A*t)v. The action of the matrix expoential on a vector
//
// @param A		The coefficient matrix for the system of ODE's
// @param v0	Initial condition vector
// @param t		Time step of the solve
//*****************************************************************************
VectorD cauchy::apply(const SparseMatrixD& A, const VectorD& v0, double t){

	// The sparse LU solver object
	Eigen::SparseLU<SparseMatrixCLD, COLAMDOrdering<int> > solver;
	
	// MPI stuff
	int myid = mpi.rank;			// Processor ID
	int numprocs = mpi.size;	// Number of processors
	int eleCount = A.rows();	// Number of species in the system

	// Number of poles
	int s = theta.rows();	
	SparseMatrixCLD At(A.rows(),A.cols());			// Scaled matrix with time
	SparseMatrixCLD tempA(A.rows(),A.cols());		// Temp variable in solution
	VectorCLD v, tempB, v0cd, myV; 
	VectorD solutionVector;								// Solution
	v0cd = v0.cast<std::complex<long double>>();
	SparseMatrixCLD ident(A.rows(),A.cols());		// Identity matrix

	myV = 0.*v0cd, v = 0*v0cd;
	ident.setIdentity();									// Builds matrix
	At = A.cast<std::complex<long double>>()*t;	// Change matrix type
	solver.analyzePattern(At);							// analyze sparsisty pattern

	// Loops over the imaginary poles. This is a linear solve over 8 lineary 
	// independent systems. The sum of all the independent solutions is w.
	for (int k = myid; k < s; k += numprocs){
		tempA = At - theta(k)*ident;
		tempB = alpha(k)*v0cd;
		// Compute the numerical factorization
		solver.factorize(tempA);

		myV = myV + solver.solve(tempB);
	}
	if (myid != 0){
		// Sends solution data to the master node 
		mpi.send(myV, eleCount, 0, MTAG1);
	}
	else {
		v = myV;
		// Receives data from the slave nodes
		for (int islave = 1; islave < numprocs; islave++) {
			myV = mpi.recv(myV, eleCount, islave, MTAG1);
			v = v + myV;
		}
	}
	v = 2.*v.real();
	v = v + alpha_0*v0cd;
	solutionVector = v.real().cast<double>();

	return solutionVector;
}

//*****************************************************************************
// Calculates exp(A*t). The matrix exponential
//
// @param A		The coefficient matrix for the system of ODE's
// @param t		Time step of the solve
//*****************************************************************************
SparseMatrixD cauchy::compute(const SparseMatrixD& A, double t){

	// MPI stuff
	int myid = mpi.rank;
	int numprocs = mpi.size;
	int eleCount = A.rows()*A.rows();

	// Number of poles
	int s = theta.rows();
	SparseMatrixCLD At(A.rows(),A.cols());
	SparseMatrixCLD tempA(A.rows(),A.cols()); 
	SparseMatrixCLD myExpA;
	SparseMatrixCLD expA;
	SparseMatrixCLD ident(A.rows(),A.cols());
	SparseMatrixLD expAReal;
	VectorD solutionVector;

	myExpA = 0*A.cast<std::complex<long double>>();
	ident.setIdentity();
	expA = 0*A.cast<std::complex<long double>>();
	At = A.cast<std::complex<long double>>()*t;

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

	return expAReal.cast<double>();
}

//*****************************************************************************
// CRAM Methods
//*****************************************************************************

//*****************************************************************************
// Initilizer for the CRAM solver
//*****************************************************************************
CRAM::CRAM(){
	MatrixCLD thetaCRAM(8,1);
	MatrixCLD alphaCRAM(8,1);

	// Defines the complex values for CRAM of order 16
	std::complex<long double> tCRAM1(-10.843917078696988026L, 19.277446167181652284L);
	std::complex<long double> tCRAM2(-5.2649713434426468895L, 16.220221473167927305L);
	std::complex<long double> tCRAM3(5.9481522689511774808L, 3.5874573620183222829L);
	std::complex<long double> tCRAM4(3.5091036084149180974L, 8.4361989858843750826L);
	std::complex<long double> tCRAM5(6.4161776990994341923L, 1.1941223933701386874L);
	std::complex<long double> tCRAM6(1.4193758971856659786L, 10.925363484496722585L);
	std::complex<long double> tCRAM7(4.9931747377179963991L, 5.9968817136039422260L);
	std::complex<long double> tCRAM8(-1.4139284624888862114L, 13.497725698892745389L);

	thetaCRAM(0,0) = tCRAM1; thetaCRAM(1,0) = tCRAM2; thetaCRAM(2,0) = tCRAM3; 
	thetaCRAM(3,0) = tCRAM4; thetaCRAM(4,0) = tCRAM5; thetaCRAM(5,0) = tCRAM6; 
	thetaCRAM(6,0) = tCRAM7; thetaCRAM(7,0) = tCRAM8;
	theta = thetaCRAM;

	// Defines the complex values for CRAM of order 16
	std::complex<long double> aCRAM1(-.0000005090152186522491565L,-.00002422001765285228797L);
	std::complex<long double> aCRAM2(.00021151742182466030907L, .0043892969647380673918L);
	std::complex<long double> aCRAM3(113.39775178483930527L, 101.9472170421585645L);
	std::complex<long double> aCRAM4(15.059585270023467528L, -5.7514052776421819979L);
	std::complex<long double> aCRAM5(-64.500878025539646595L, -224.59440762652096056L);
	std::complex<long double> aCRAM6(-1.4793007113557999718L, 1.7686588323782937906L);
	std::complex<long double> aCRAM7(-62.518392463207918892L, -11.19039109428322848L);
	std::complex<long double> aCRAM8(.041023136835410021273L, -.15743466173455468191L);

	alphaCRAM(0,0) = aCRAM1; alphaCRAM(1,0) = aCRAM2; alphaCRAM(2,0) = aCRAM3; 
	alphaCRAM(3,0) = aCRAM4; alphaCRAM(4,0) = aCRAM5; alphaCRAM(5,0) = aCRAM6; 
	alphaCRAM(6,0) = aCRAM7; alphaCRAM(7,0) = aCRAM8; 
	alpha_0 = 2.1248537104952237488e-16L;
	alpha = alphaCRAM;
}

//*****************************************************************************
// Parabolic Methods
//*****************************************************************************

//*****************************************************************************
// Initilizer for the parabolic solver
//*****************************************************************************
parabolic::parabolic(){
	// Gets the coefficients for the contour integral
	MatrixCLD coeffs = parabolicContourCoeffs(order);
	theta = coeffs.col(0);
	alpha = coeffs.col(1);
}

//*****************************************************************************
// Calcualtes the quadrature points for the parabolic contour
//
// @param order	Order of the quadrature approx
//*****************************************************************************
MatrixCLD parabolic::parabolicContourCoeffs(int N){
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
// Hyperbolic Methods
//*****************************************************************************

//*****************************************************************************
// Initilizer for the hyperbolic solver
//*****************************************************************************
hyperbolic::hyperbolic(){
	// Gets the coefficients for the contour integral
	MatrixCLD coeffs = hyperbolicContourCoeffs(order);
	theta = coeffs.col(0);
	alpha = coeffs.col(1);
}

//*****************************************************************************
// Calcualtes the quadrature points for the hyperbolic contour
//
// @param order	Order of the quadrature approx
//*****************************************************************************
MatrixCLD hyperbolic::hyperbolicContourCoeffs(int N){
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
