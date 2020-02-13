#include "matrixExponential.h"
#include <unsupported/Eigen/MatrixFunctions>
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
	else if (type == "pade-method1"){
		solver = new method1;
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
// Methods for Pade Class
//
// The pade functions are taken from the eigen3 unsuported matrix exponential
// code. I didn't like the linear solver they used and the solver they used
// is not avaliable for sparsematrices so idk how it work when you pass in one.
//*****************************************************************************

//*****************************************************************************
// Pade approximation for order (3,3)
//
// @param A		Sparse matrix
// @param A2	Sparse matrix A^2
//*****************************************************************************
void pade::pade3(const SparseMatrixD& A, const SparseMatrixD& A2, 
	SparseMatrixD& U, SparseMatrixD& V){
	// Pade coefficients
	long double b[] = {120.L, 60.L, 12.L, 1.L};
	SparseMatrixD ident(A.rows(), A.cols()), temp;

	ident.setIdentity();
	temp = b[3]*A2 + b[1]*ident;
	U = A*temp;
	V = b[2]*A2 + b[0]*ident;

}

//*****************************************************************************
// Pade approximation for order (5,5)
//
// @param A		Sparse matrix
// @param A2	Sparse matrix A^2
// @param A4	Sparse matrix A^4
//*****************************************************************************
void pade::pade5(const SparseMatrixD& A, const SparseMatrixD& A2, 
	const SparseMatrixD& A4, SparseMatrixD& U, SparseMatrixD& V){
	// Pade coefficients
	long double b[] = {30240.L, 15120.L, 3360.L, 420.L, 30.L, 1.L};
	SparseMatrixD ident(A.rows(), A.cols()), temp;

	ident.setIdentity();
	temp = b[5]*A4 + b[3]*A2 + b[1]*ident;
	U = A*temp;
	V = b[4]*A4 + b[2]*A2 + b[0]*ident;
	
}

//*****************************************************************************
// Pade approximation for order (7,7)
//
// @param A		Sparse matrix
// @param A2	Sparse matrix A^2
// @param A4	Sparse matrix A^4
// @param A6	Sparse matrix A^6
//*****************************************************************************
void pade::pade7(const SparseMatrixD& A, const SparseMatrixD& A2, 
	const SparseMatrixD& A4, const SparseMatrixD& A6, SparseMatrixD& U, 
	SparseMatrixD& V){
	// Pade coefficients
	long double b[] = {17297280.L, 8648640.L, 1995840.L, 277200.L, 25200.L, 
							 1512.L, 56.L, 1.L};
	SparseMatrixD ident(A.rows(), A.cols()), temp;

	ident.setIdentity();
	temp = b[7] * A6 + b[5]*A4 + b[3]*A2 + b[1]*ident;
	U = A*temp;
	V = b[4]*A4 + b[2]*A2 + b[0]*ident;
	
}

//*****************************************************************************
// Pade approximation for order (9,9)
//
// @param A		Sparse matrix
// @param A2	Sparse matrix A^2
// @param A4	Sparse matrix A^4
// @param A6	Sparse matrix A^6
// @param A8	Sparse matrix A^8
//*****************************************************************************
void pade::pade9(const SparseMatrixD& A, const SparseMatrixD& A2, 
	const SparseMatrixD& A4, const SparseMatrixD& A6, const SparseMatrixD& A8,
	SparseMatrixD& U, SparseMatrixD& V){
	// Pade coefficients
	long double b[] = {17643225600.L, 8821612800.L, 2075673600.L, 302702400.L, 
						    30270240.L, 2162160.L, 110880.L, 3960.L, 90.L, 1.L}; 
	SparseMatrixD ident(A.rows(), A.cols()), temp;
	ident.setIdentity();
	temp = b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident;
	U = A*temp;
	V = b[8]*A8 + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident;

}

//*****************************************************************************
// Pade approximation for order (13,13)
//
// @param A		Sparse matrix
// @param A2	Sparse matrix A^2
// @param A4	Sparse matrix A^4
// @param A6	Sparse matrix A^6
//*****************************************************************************
void pade::pade13(const SparseMatrixD& A, const SparseMatrixD& A2, 
	const SparseMatrixD& A4, const SparseMatrixD& A6, SparseMatrixD& U, 
	SparseMatrixD& V){
	// Pade coefficients
	long double b[] = {64764752532480000.L, 32382376266240000.L, 
							 7771770303897600.L, 1187353796428800.L, 129060195264000.L, 
							 10559470521600.L, 670442572800.L, 33522128640.L, 
							 1323241920.L, 40840800.L, 960960.L, 16380.L, 182.L, 1.L};
	SparseMatrixD ident(A.rows(), A.cols()), temp;
	ident.setIdentity();
	V = b[13]*A6 + b[11]*A4 + b[9]*A2;
	temp = A6 * V;
	temp += b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident;
	U = A*temp;
	temp = b[12]*A6 + b[10]*A4 + b[8]*A2;
	V = A6*temp;
	V += b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident;
}
//*****************************************************************************
// Calculates exp(A*t)v. The action of the matrix expoential on a vector
//
// @param A		The coefficient matrix for the system of ODE's
// @param v0	Initial condition vector
// @param t		Time step of the solve
//*****************************************************************************
VectorD pade::apply(const SparseMatrixD& A, const VectorD& v0, double t){
	SparseMatrixD matExp = compute(A, t);
	return matExp*v0;
}

//*****************************************************************************
// Calculates exp(A*t). The the matrix expoential
//
// @param A		The coefficient matrix for the system of ODE's
// @param t		Time step of the solve
//*****************************************************************************
SparseMatrixD pade::compute(const SparseMatrixD& A, double t){
	// The sparse LU solver object
	Eigen::SparseLU<SparseMatrixD, COLAMDOrdering<int> > solver;
	SparseMatrixD U, V, At, denominator, numerator, R;
	int alpha;
	// Calculate At matrix
	At = A*t;
	// Runs the pade algorithm to find matrices U and V
	run(At, U, V, alpha);	
	// Builds numerator and denominator
	denominator = -U + V;
	numerator = U + V;
	// Analyze sparcisity patern
	solver.analyzePattern(denominator);
	// Compute LU decomp
	solver.factorize(denominator);

	// Solve for the matrix exponential
	R = solver.solve(numerator);
	// unscale the matrix if needed
	for (int k=0; k<alpha; k++){
		R = R*R;
	}
	return R;	
}

//*****************************************************************************
// Run method for pade class method1
//
// @param A			Sparse matrix
// @param U			U matrix of pade solver
// @param V			V matrix of pade solver
// @param alpha	Number of times to square the matrix
//*****************************************************************************
void method1::run(const SparseMatrixD& A, SparseMatrixD& U, SparseMatrixD& V,
	int& alpha){
	const double l1Norm = (A.cwiseAbs()*VectorD::Ones(A.cols())).maxCoeff();
	SparseMatrixD Ascaled;
	double maxnorm, scale;
	const SparseMatrixD A2 = A*A;
	alpha = 0;

	if (l1Norm < 1.495585217958292e-002){
		pade3(A, A2, U, V);
	}
	else if (l1Norm < 2.539398330063230e-001){
		const SparseMatrixD A4 = A2*A2;
		pade5(A, A2, A4, U, V);
	}
	else if (l1Norm < 9.504178996162932e-001){
		const SparseMatrixD A4 = A2*A2;
		const SparseMatrixD A6 = A4*A2;
		pade7(A, A2, A4, A6, U, V);
	}
	else if (l1Norm < 2.097847961257068e+000){
		const SparseMatrixD A4 = A2*A2;
		const SparseMatrixD A6 = A4*A2;
		const SparseMatrixD A8 = A6*A2;
		pade9(A, A2, A4, A6, A8, U, V);
	}
	else{
		maxnorm = 5.371920351148152;
		// Calculate the number of squarings
		std::frexp(l1Norm/maxnorm, &alpha);
		if (alpha < 0){alpha = 0;};
		// Get the scale
		scale = std::pow(2,alpha); 
		// Scale the matrix
		Ascaled = A/scale;
		// Calculate the squared matrices based on the new scaled matrix
		const SparseMatrixD A2scaled = Ascaled*Ascaled;
		const SparseMatrixD A4scaled = A2scaled*A2scaled;
		const SparseMatrixD A6scaled = A4scaled*A2scaled;
		pade13(Ascaled, A2scaled, A4scaled, A6scaled, U, V);
	}
}


//*****************************************************************************
// Runs the algorithm for method 2
//
// @param A			Sparse matrix
// @param U			U matrix of pade solver
// @param V			V matrix of pade solver
// @param alpha	Number of times to square the matrix
//*****************************************************************************
void method2::run(const SparseMatrixD& A, SparseMatrixD& U, SparseMatrixD& V,
	int& alpha){
	double l1Norm = (A.cwiseAbs()*VectorD::Ones(A.cols())).maxCoeff();
}

//*****************************************************************************
// Normest
// Produces the l1norm of A*B
//
// @param A		Sparse matrix
// @param B		Sparse matrix
//*****************************************************************************
double method2::normest(const SparseMatrixD& A, const SparseMatrixD& B){
	const SparseMatrixD& C = A*B;
	double l1Norm = (C.cwiseAbs()*VectorD::Ones(A.cols())).maxCoeff();
	return l1Norm;
}

//*****************************************************************************
// Normest
// Produces the l1norm of A^m
//
// @param A		Sparse matrix
// @param m		integer, power of the matrix
//*****************************************************************************
double method2::normest(const SparseMatrixD& A, const int m){
	SparseMatrixD C(A.rows(), A.cols());
	double l1Norm;

	// Rises the matrix to power m
	for (int i; i<m; i++){
		C = C*A;
	}
	l1Norm = (C.cwiseAbs()*VectorD::Ones(A.cols())).maxCoeff();
	return l1Norm;
}
//*****************************************************************************
// Ell
// Returns the integer max((log2(alpha/u)/(2m)), 0), where 
// alpha = = |c_(2m+1)|normest(|A|, 2m + 1)/l1norm(A1)
//
// @param A		Sparse matrix
// @param m		integer, power of the matrix
//*****************************************************************************
int method2::ell(const SparseMatrixD& A, const int m){
	double c, u, A1NormMatrixPower, A1Norm, alpha, log2AlphaOverU;
	int p, value;
	p = 2*m + 1;

	// Leading coeffcient
	c = std::abs((double)1./(binomialCoeff(2*p, p)*factorial(2*p+1))); 
	// unit round off IEE double
	u = std::pow(2,-53); 
	// l1 norm of matrix power
	A1NormMatrixPower = normest(A.cwiseAbs(), p);
	// l1 norm of matrix
	A1Norm = (A.cwiseAbs()*VectorD::Ones(A.cols())).maxCoeff();	

	alpha = c*A1NormMatrixPower/A1Norm;
	log2AlphaOverU = std::log2(alpha/u);
	value = (int) std::ceil(log2AlphaOverU/(2.*m));
	return std::max(value, 0);
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
	const int myid = mpi.rank;			// Processor ID
	const int numprocs = mpi.size;	// Number of processors
	const int eleCount = A.rows();	// Number of species in the system

	// Number of poles
	
	const int s = theta.rows();	
	SparseMatrixCLD At(A.rows(),A.cols());		// Scaled matrix with time
	SparseMatrixCLD tempA(A.rows(),A.cols());			// Temp variable in solution
	VectorCLD v, tempB, v0cd, myV; 
	VectorD solutionVector;									// Solution
	v0cd = v0.cast<std::complex<long double>>();
	SparseMatrixCLD ident(A.rows(),A.cols());	// Identity matrix

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
	const int myid = mpi.rank;
	const int numprocs = mpi.size;
	const int eleCount = A.rows()*A.rows();

	// Number of poles
	const int s = theta.rows();
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
