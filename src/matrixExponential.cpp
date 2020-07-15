#include "matrixExponential.h"
//*****************************************************************************
// Methods for MatrixExponential factory class
//*****************************************************************************
matrixExponential *matrixExponentialFactory::getExpSolver(std::string type,
	bool krylovBool, int krylovDim){
	matrixExponential *solver = nullptr;

	if (type == "CRAM"){
		solver = new CRAM(krylovBool, krylovDim);
		return solver;
	}
	else if (type == "parabolic"){
		solver = new parabolic(krylovBool, krylovDim);
		return solver;
	}
	else if (type == "hyperbolic"){
		solver = new hyperbolic(krylovBool, krylovDim);
		return solver;
	}
	else if (type == "pade-method1"){
		solver = new method1(krylovBool, krylovDim);
		return solver;
	}
	else if (type == "pade-method2"){
		solver = new method2(krylovBool, krylovDim);
		return solver;
	}
	else if (type == "LPAM"){
		solver = new LPAM(krylovBool, krylovDim);
		return solver;
	}
	else if (type == "taylor"){
		solver = new taylor(krylovBool, krylovDim);
		return solver;
	}
	else {
		std::string errorMessage =
			" You have selected a matrix exponential solver\n"
			" that is not in libowski. Avaliable solvers are\n\n"
			" CRAM\n"
			" parabolic\n"
			" hyperbolic\n"
			" pade-method1\n"
			" pade-method2\n"
			" LPAM\n"
			" Taylor\n";
		libowskiException::runtimeError(errorMessage);
		return solver;
	}
}


//*****************************************************************************
// Default constructor for matrixExponential class
//
// @param KrylovFlag		set to true if you want to use the Krylov subspace 
//								in the calculation. WARNING THIS CAN ONLY BE USED WITH
//								THE APPLY FUNCTION
//	@param subspaceDim	The dimension of the krylov subspace
//*****************************************************************************
matrixExponential::matrixExponential(bool KrylovFlag, int subspaceDim){

	useKrylovSubspace = KrylovFlag;
	krylovSubspaceDim = subspaceDim;
}

//*****************************************************************************
// Methods for Taylor series Class
//
//*****************************************************************************
taylor::taylor(bool krylovBool, int krylovDim):matrixExponential(krylovBool,
	krylovDim){
	std::string thetaFname = getDataPath() + "theta.txt";
   std::ifstream inFile;
   double x;
   int counter = 0;
	name = "taylor";

   inFile.open(thetaFname);
   if (!inFile){
		std::string errorMessage =
			" Unable to find the Taylor solver theta file\n";
		libowskiException::runtimeError(errorMessage);
   }   
   while (inFile >> x) {
      theta(counter) = x;
      counter ++; 
   }   	
}

//*****************************************************************************
// Selects the Taylor series degree for the approximation
//
// @param A				Sparse transition matrix
// @param b				Initional condition matrix
// @param M				Returned matrix M
// @param m_max		Max order of Taylor series expansion
// @param p_max		Parameter from paper
// @param shift		Bool to shift the norm of the matrix
// @param forceEstm	Bool for the hoice of how to choose m
//*****************************************************************************
void taylor::parameters(const SparseMatrixD& A, const VectorD& b, MatrixD& M,
	int m_max, int p_max, bool shift, bool forceEstm){
	int n = A.cols();
	double normA, tempLim, c, mu;
	MatrixD eta, alpha;
	SparseMatrixD ident(A.rows(), A.cols()), Ai;
	tempLim = 4.*theta(m_max)*(double)p_max*((double)p_max+3.)/
		((double)m_max*(double)b.rows());

	// A matrix that is internal to the fucntion
	Ai = A;
	ident.setIdentity();

	if (shift){
		mu = A.diagonal().sum()/((double)n);
		Ai = Ai - mu*ident;
	}
	if (not forceEstm){ normA = l1norm(Ai);};
	if (not forceEstm and normA <= tempLim){
		c = normA;
		alpha = c*MatrixD::Ones(p_max-1,1);
	}
	else{
		eta = MatrixD::Zero(p_max,1);
		alpha = MatrixD::Zero(p_max-1,1);
		for (int p = 1; p < p_max+1; p++){
			c = normAm(Ai, p+1);
			c = std::pow(c, 1./((double)p+1.));
			eta(p-1) = c;
		}

		for (int p = 1; p < p_max; p++){
			alpha(p-1) = std::max(eta(p-1), eta(p));
		}
	}
	M = MatrixD::Zero(m_max, p_max-1);
	for (int p = 2; p < p_max+1; p++){
		for (int m = p*(p-1)-1; m < m_max+1; m++){
			M(m-1, p-2) = alpha(p-2)/theta(m-1);
		}
	}
}
//*****************************************************************************
// Private version
// Calculates exp(A*t)v. The action of the matrix expoential on a vector
//*****************************************************************************
VectorD taylor::expmv(const SparseMatrixD& A, const double t, const VectorD& v0,
	MatrixD& M, bool shift, bool fullTerm){
	int n = A.cols(), m, p, m_max;
	long int cost = 1e10L, s = 1L;
	double tol = std::pow(2.,-53.), tt, mu, eta, c1, c2;
	SparseMatrixD ident(A.rows(), A.cols()), Ai = A;
	VectorLI diag, minCols;
	MatrixLI C, U;
	VectorD f, b = v0;

	ident.setIdentity();
	if (shift){
		mu = Ai.diagonal().sum()/((double)n);
		Ai = Ai - mu*ident;
	}

	if (M.size() == 0){
		tt = 1.;
		parameters(Ai*t, b, M);
	}
	else{
		tt = t;
	}

	if (t == 0){
		m = 0;
	}
	else{
		m_max = M.rows();
		p = M.cols();
		diag = VectorLI::Zero(m_max);
		U = MatrixLI::Zero(m_max, m_max);
		for (int i=1; i<m_max+1; i++){ diag(i-1) = i;};
		U.diagonal() = diag;
		C = MatrixLI((std::abs(tt)*M).array().ceil().cast<long int>()).transpose() * U;
		minCols = C.colwise().minCoeff();
		for (int i = 1; i < minCols.rows(); i++){ 
			if (std::labs(minCols(i)) < cost and std::labs(minCols(i)) > 0){
				cost = std::labs(minCols(i));
				m = i;
			}
		}
		m = m + 1;
		s = std::max(cost/m, 1L);
	}
	if (shift){
		eta = std::exp(t*mu/(double)s);
	}
	else{
		eta = 1.;
	}
	f = b;
	for (int i = 1; i < s+1; i++){
		c1 = b.lpNorm<Infinity>();
		for (int k = 1; k	< m+1; k++){
			b = (t/((double)s*(double)k))*(Ai*b);
			f = f + b;
			c2 = b.lpNorm<Infinity>();
			if (not fullTerm){
				if (c1 + c2 <= tol*f.lpNorm<Infinity>()){
					break;
				}
			}
		}
		f = eta*f;
		b = f;
	}
	return f;		 
}

//*****************************************************************************
// Calculates exp(A*t)v. The action of the matrix expoential on a vector
//
// @param A		The coefficient matrix for the system of ODE's
// @param v0	Initial condition vector
// @param t		Time step of the solve
//*****************************************************************************
VectorD taylor::apply(const SparseMatrixD& A, const VectorD& v0, double t){
	MatrixD M;
	VectorD sol;
	sol = expmv(A, t, v0, M);
	return sol;
}

//*****************************************************************************
// Calculates exp(A*t). The the matrix expoential
//
// @param A		The coefficient matrix for the system of ODE's
// @param t		Time step of the solve
//*****************************************************************************
SparseMatrixD taylor::compute(const SparseMatrixD& A, double t){
	std::string errorMessage =
		" The Taylor series method cannot be used with compute\n";
	libowskiException::runtimeError(errorMessage);
	return A;
}

//*****************************************************************************
// Methods for Pade Class
//
// The pade functions are taken from the eigen3 unsuported matrix exponential
// code. I didn't like the linear solver they used and the solver they used
// is not avaliable for sparsematrices so idk how it work when you pass in one.
//*****************************************************************************

//*****************************************************************************
// Constructer for pade class
//*****************************************************************************
pade::pade(bool krylovBool, int krylovDim):matrixExponential(krylovBool, 
	krylovDim){};

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
	V = b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident;
	
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
	SparseMatrixD matExp;
	// Generates the krylov subspace if needed
	if (useKrylovSubspace){
		// define matrices and variables
		MatrixD Q;
		SparseMatrixD H;
		SparseMatrixD At = A*t;
		double beta = v0.norm();
		// Builds the krylov subspace
		arnoldi(At, v0, krylovSubspaceDim, Q, H);
		// Computes matrix exponential of krylov subsapce
		matExp = compute(H, 1.0);	
		// Computes the result
		return beta*(Q*matExp*VectorD::Unit(krylovSubspaceDim,0));
	}
	// Computes it without krylov subspace
	matExp = compute(A, t);
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
	SparseMatrixD U, V, H, At, denominator, numerator, R;
	MatrixD Q;
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
// Constructor for method1
//*****************************************************************************
method1::method1(bool krylovBool, int krylovDim):pade(krylovBool, krylovDim){
	name = "pade-method1";
};

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
	const double A1norm = l1norm(A);
	SparseMatrixD Ascaled;
	double maxnorm, scale;
	const SparseMatrixD A2 = A*A;
	alpha = 0;

	if (A1norm < 1.495585217958292e-002){
		pade3(A, A2, U, V);
	}
	else if (A1norm < 2.539398330063230e-001){
		const SparseMatrixD A4 = A2*A2;
		pade5(A, A2, A4, U, V);
	}
	else if (A1norm < 9.504178996162932e-001){
		const SparseMatrixD A4 = A2*A2;
		const SparseMatrixD A6 = A4*A2;
		pade7(A, A2, A4, A6, U, V);
	}
	else if (A1norm < 2.097847961257068e+000){
		const SparseMatrixD A4 = A2*A2;
		const SparseMatrixD A6 = A4*A2;
		const SparseMatrixD A8 = A6*A2;
		pade9(A, A2, A4, A6, A8, U, V);
	}
	else{
		maxnorm = 5.371920351148152;
		// Calculate the number of squarings
		std::frexp(A1norm/maxnorm, &alpha);
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
// Constructor for method2
//*****************************************************************************
method2::method2(bool krylovBool, int krylovDim):pade(krylovBool, krylovDim){
	name = "pade-method2";
};

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
	double d4, d6, d8, d10, eta1, eta2, eta3, eta4, eta5, log2EtaOverTheta;
	int value;
	SparseMatrixD B, B2, B4, B6;
	const double A1norm = l1norm(A);
	const double theta13 = 4.25;
	alpha = 0;

	// try pade3
	const SparseMatrixD A2 = A*A;
	d6 = std::pow(normAm(A2,3), 1./6.);
	eta1 = std::max(std::pow(normAm(A2,2), 1./4.), d6);
	if (eta1 < 1.495585217958292e-002 and ell(A,3) == 0){
		// Calculate U and V using pade3
		pade3(A, A2, U, V);
		// exit the function
		return;
	}

	// try pade5
	const SparseMatrixD A4 = A2*A2;
	d4 = std::pow(l1norm(A4), 1./4.);
	eta2 = std::max(d4, d6);
	if (eta2 < 2.539398330063230e-001 and ell(A,5) == 0){
		// Calculate U and V using pade5
		pade5(A, A2, A4, U, V);
		// exit the function
		return;

	}
	// try pade 7 and 9
	const SparseMatrixD A6 = A4*A2;
	d6 = std::pow(l1norm(A6), 1./6.);
	d8 = std::pow(normAm(A4,2),1./8.);
	eta3 = std::max(d6, d8);
	if (eta3 < 9.504178996162932e-001 and ell(A,7) == 0){
		// Calculate U and V using pade7
		pade7(A, A2, A4, A6, U, V);
		// exit the function
		return;
	}
	if (eta3 < 2.097847961257068e+000 and ell(A,9) == 0){
		const SparseMatrixD A8 = A6*A2;
		// Calculate U and V using pade9
		pade9(A, A2, A4, A6, A8, U, V);
		// exit the function
		return;
	}

	// Do pade 13
	d10 = std::pow(normAm(A4, A6), 1./10.);
	eta4 = std::max(d8, d10);
	eta5 = std::min(eta3, eta4);
	log2EtaOverTheta = std::log2(eta5/theta13);
	value = (int) std::ceil(log2EtaOverTheta);
	// find number of squarings 
	alpha = std::max(value, 0);
	alpha = alpha + ell(A*std::pow(2,-alpha), 13);
	// Scale down the matrices
	B = A*std::pow(2,-alpha);
	B2 = A2*std::pow(2,-2*alpha);
	B4 = A4*std::pow(2,-4*alpha);
	B6 = A6*std::pow(2,-6*alpha);
	// Solve for U and V
	pade13(B, B2, B4, B6, U, V);
	// exit function
	return;

}

//*****************************************************************************
// Ell
// Returns the integer max((log2(alpha/u)/(2m)), 0), where 
// alpha = = |c_(2m+1)|normAm(|A|, 2m + 1)/l1norm(A1)
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
	A1NormMatrixPower = normAm(A.cwiseAbs(), p);
	// l1 norm of matrix
	A1Norm = l1norm(A);

	alpha = c*A1NormMatrixPower/A1Norm;
	log2AlphaOverU = std::log2(alpha/u);
	value = (int) std::ceil(log2AlphaOverU/(2.*m));
	return std::max(value, 0);
}

//*****************************************************************************
// Methods for Cauchy Class
//*****************************************************************************

//*****************************************************************************
// Cauchy constructor
//*****************************************************************************
cauchy::cauchy(bool krylovBool, int krylovDim):matrixExponential(krylovBool, 
	krylovDim){};

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

	// Number of poles
	const int s = theta.rows();	
	SparseMatrixCLD At(A.rows(),A.cols());		// Scaled matrix with time
	SparseMatrixCLD tempA(A.rows(),A.cols());			// Temp variable in solution
	VectorCLD v, tempB, v0cd, myV; 
	VectorD solutionVector;									// Solution
	SparseMatrixCLD ident(A.rows(),A.cols());	// Identity matrix

	ident.setIdentity();									// Builds matrix
	v0cd = v0.cast<std::complex<long double>>(); // Change the vector type
	myV = 0.*v0cd, v = 0.*v0cd;
	At = A.cast<std::complex<long double>>()*t;	// Change matrix type
	solver.analyzePattern(At);							// analyze sparsisty pattern

	// Sends the correct matrix and vector to each processor if the processor
	// isn't the main node
	if (mpi.rank == 0){ 
		// Sends data from the worker nodes
		for (int worker = 1; worker < numprocs; worker++) {
			mpi.send(At, worker, 100);
			mpi.send(v0cd, worker, 200);
		}
	}
	else{
		// Receives the data from the main node
		mpi.recv(At, 0, 100);
		mpi.recv(v0cd, 0, 200);
	}

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
		// Sends solution data to the main node 
		mpi.send(myV, 0, 300);
	}
	else {
		v = myV;
		// Receives data from the worker nodes
		for (int worker = 1; worker < numprocs; worker++) {
			mpi.recv(myV, worker, 300);
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

	// If the Krylov subspace flag is true kill the program
	if (useKrylovSubspace){
		std::string errorMessage = " Krylov subspace method cannot \n"
			" be used with the compute method\n";
		libowskiException::runtimeError(errorMessage);
	}

	// MPI stuff
	const int myid = mpi.rank;
	const int numprocs = mpi.size;

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

	// Sends the correct matrix and vector to each processor if the processor
	// isn't the main node
	if (mpi.rank == 0){ 
		// Sends data from the worker nodes
		for (int worker = 1; worker < numprocs; worker++) {
			mpi.send(At, worker, 100);
		}
	}
	else{
		// Receives the data from the main node
		mpi.recv(At, 0, 100);
	}

	// Loops over the imaginary poles to compute the matrix exponential
	for (int k = myid; k < s; k += numprocs){
		tempA = At - theta(k)*ident;

		// Multipy A-1*alpha(K) because the matrix is 1-1 and only left invertible 
		// because it has linearly independent columns but not rows. 
		myExpA = myExpA + MoorePenroseInv(tempA)*alpha(k);
	}
	if (myid != 0){
		// Sends solution data to the main node 
		mpi.send(myExpA, 0, 200);
	}
	else {
		expA = myExpA;
		// Receives data from the worker nodes
		for (int worker = 1; worker < numprocs; worker++) {
			mpi.recv(myExpA, worker, 200);
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
CRAM::CRAM(bool krylovBool, int krylovDim):cauchy(krylovBool, krylovDim
	){

	name = "CRAM";
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
parabolic::parabolic(bool krylovBool, int krylovDim):cauchy(krylovBool, 
	krylovDim){

	name = "parabolic";
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
hyperbolic::hyperbolic(bool krylovBool, int krylovDim):cauchy(krylovBool, 
	krylovDim){

	name = "hyperbolic";
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

//*****************************************************************************
// LPAM Methods
//*****************************************************************************

//*****************************************************************************
// Initilizer for the LPAM solver
//*****************************************************************************
LPAM::LPAM(bool krylovBool, int krylovDim):matrixExponential(krylovBool, 
	krylovDim){};

//*****************************************************************************
// Calculates the leading coefficient for the Laguerre Polynomial
//
// @param tau	Scaling coefficient for Leguerre Polynomial
// @param a		Some other factor thing that is used for Leguerre Polynomials
// @param n		Coefficient of expansion
//*****************************************************************************
double LPAM::laguerreCoefficient(double tau, double a, double n){
	double temp1, temp2;	

	if (n == 0){
		return 1.;
	}
	else if (n == 1){
		return 1. + a - tau;
	}
	else{
		temp1 = (2.*(n - 1.) + 1. + a - tau)*laguerreCoefficient(tau, a, n-1);
		temp2 = (n - 1. + a)*laguerreCoefficient(tau, a, n-2);
		return 1./(n)*(temp1 - temp2);
	}

}

//*****************************************************************************
// Calculates the Laguerre Polynomial vector 
//
// @param tau		Scaling coefficient for Leguerre Polynomial
// @param a			Some other factor thing that is used for Leguerre Polynomials
// @param k			The current sum iteration
// @param n0		Initial condition
// @param AScaled	The scaled matrix
//*****************************************************************************
VectorD LPAM::laguerrePolynomial(double tau, int a, int k, const VectorD& n0, 
	const SparseMatrixD& AScaled){
	double lc = std::pow(-1., k);
	VectorD solveProduct, vect;
	SparseMatrixD ident(AScaled.rows(),AScaled.cols());
	ident.setIdentity();

	vect = 0.*n0;
	solveProduct = solver.solve(n0);

	// Build the matrix inverse part of the Polynomial
	for (int n = 0; n < a+k; n++){
		solveProduct = solver.solve(solveProduct);
	}

	if (k == 0){ return ident * solveProduct;};

	vect = AScaled * solveProduct;

	for (int n = 0; n < k-1; n++){
		vect = AScaled * vect;
	}

	return lc*vect;
}

//*****************************************************************************
// Calculates the Laguerre Polynomial matrix
//
// @param tau		Scaling coefficient for Leguerre Polynomial
// @param a			Some other factor thing that is used for Leguerre Polynomials
// @param k			The current sum iteration
// @param AScaled	The scaled matrix
//*****************************************************************************
SparseMatrixD LPAM::laguerrePolynomial(double tau, int a, int k, const 
	SparseMatrixD& AScaled){
	double lc = std::pow(-1., k);
	SparseMatrixD solveProduct, matrix;
	SparseMatrixD ident(AScaled.rows(),AScaled.cols());
	ident.setIdentity();

	matrix = 0.*AScaled;
	solveProduct = matInverse;
	if (k == 0){ return ident * solveProduct;};

	// Build the matrix inverse part of the Polynomial
	for (int n = 0; n < a+k; n++){
		solveProduct = solveProduct * matInverse;
	}

	matrix = AScaled * solveProduct;

	for (int n = 0; n < k-1; n++){
		matrix = AScaled * matrix;
	}

	return lc*matrix;
}

//*****************************************************************************
// Calculates exp(A*t)v. The action of the matrix expoential on a vector
//
// @param A		The coefficient matrix for the system of ODE's
// @param v0	Initial condition vector
// @param t		Time step of the solve
//*****************************************************************************
VectorD LPAM::apply(const SparseMatrixD& A, const VectorD& v0, double t){
	SparseMatrixD AScaled = A*t/tau;
	SparseMatrixD ident(AScaled.rows(),AScaled.cols());
	SparseMatrixD tempA;
	VectorD solution = 0.*v0, lp;
	double lc;
	ident.setIdentity();
	tempA = ident - AScaled;

	// Compute LU decomp
	solver.compute(tempA);

	for (int n = 0; n < k; n++){
		lc = laguerreCoefficient(tau, a, n);
		lp = laguerrePolynomial(tau, (int)a, n, v0, AScaled);
		solution = solution + lc*lp;
	}
	return solution;
}

//*****************************************************************************
// Calculates exp(A*t). The matrix expoential
//
// @param A		The coefficient matrix for the system of ODE's
// @param t		Time step of the solve
//*****************************************************************************
SparseMatrixD LPAM::compute(const SparseMatrixD& A, double t){
	SparseMatrixD AScaled = A*t/tau;
	SparseMatrixD ident(AScaled.rows(),AScaled.cols());
	SparseMatrixD tempA = ident - AScaled;
	SparseMatrixD matExp = 0.*A, lp;
	double lc;
	ident.setIdentity();

	// Compute the inverse
	matInverse = MoorePenroseInv(tempA);

	for (int n = 0; n < k; n++){
		lc = laguerreCoefficient(tau, a, n);
		lp = laguerrePolynomial(tau, (int)a, n, AScaled);
		matExp = matExp + lc*lp;
	}
	return matExp;
}
