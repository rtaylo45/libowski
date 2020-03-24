//*****************************************************************************
// Author: Zack Taylor
// 
// Classical numerical integrator class that implements Runge-Kutta
// methods. These methods solve the general form,
//*****************************************************************************
#include "ODEintegrator.h"

//*****************************************************************************
// Default constructor for ode integrator class
//
// @param solverName		Name of the solver
//*****************************************************************************
ODEintegrator::ODEintegrator(std::string solverName){
	name = solverName;
	isInit = true;
}

//*****************************************************************************
// Constructor for Runge-Kutta integrator class
//
// @param solverName		Name of the solver
// @param aCoeffs			a coefficients for Runge-Kutta
// @param bCoeffs			b coefficients for Runge-Kutta
//*****************************************************************************
rungeKuttaIntegrator::rungeKuttaIntegrator(std::string solverName, ArrayD aCoeffs,
	ArrayD bCoeffs):ODEintegrator(solverName){
	order = aCoeffs.cols();
	a = aCoeffs;
	b = bCoeffs;
	assert(order > 0);
}
//*****************************************************************************
// Constructor for explicit Runge-Kutta integrator class
//
// @param solverName		Name of the solver
// @param aCoeffs			a coefficients for Runge-Kutta
// @param bCoeffs			b coefficients for Runge-Kutta
//*****************************************************************************
explicitRKIntegrator::explicitRKIntegrator(std::string solverName, ArrayD aCoeffs,
	ArrayD bCoeffs):rungeKuttaIntegrator(solverName, aCoeffs, bCoeffs){};

//*****************************************************************************
// Methods for the explicit integrator
//
// @param A		ODE matrix L that is used in the function
// @param yn	Current species vector solution
// @param dt	Time step to integrate over
//*****************************************************************************
VectorD explicitRKIntegrator::integrate(const SparseMatrixD& A, const VectorD& yn,
	double dt){
	VectorD ynPlus1 = 0*yn, ySum = 0*yn;

	// Loops over the order	
	for (int i=1; i <= order; i++){
		ySum += b(i-1)*kn(i, A, yn, dt);	
	}
	ynPlus1 = yn + dt*ySum;
	return ynPlus1;
}

//*****************************************************************************
// Computes the k function used in Rung-Kutta
//
// @param order	Order of the k function
// @param A			ODE matrix L that is used in the function
// @param yn		Current species vector solution
// @param dt		Time step to integrate over
//*****************************************************************************
VectorD explicitRKIntegrator::kn(int order, const SparseMatrixD& A, const
	VectorD& yn, double dt){
	VectorD sum = 0*yn, tempy = 0*yn;

	// Sums over kn's to build the current kn	
	for (int j=1; j < order; j++){
		sum += a(order-1,j-1)*kn(j, A, yn, dt);
	}
	tempy = yn + dt*sum;
	return A*tempy;
}

//*****************************************************************************
// Cleans the integrator
//*****************************************************************************
void explicitRKIntegrator::clean(){};
//*****************************************************************************
// Constructor for BDF integrator
//
// @param solverName		Name of the solver
// @param BDForder		Order of the solver
//*****************************************************************************
BDFIntegrator::BDFIntegrator(std::string solverName, int BDForder):ODEintegrator(
	solverName){
	order = BDForder;
	assert(order > 0);
	MatrixD aCoeffs(6,6);
	VectorD bCoeffs(6);
	aCoeffs << 1., 0., 0., 0., 0., 0.,
		4./3., -1./3., 0., 0., 0., 0.,
		18./11., -9./11., 2./11., 0., 0., 0.,
		48./25., -36./25., 16./25., -3./25., 0., 0.,
		300./137., -300./137., 200./137., -75./137., 12./137., 0.,
		360./147., -450./147., 400./147., -225./147., 72./147., -10./147.;

	bCoeffs << 1., 2./3., 6./11., 12./25., 60./137., 60./147;
	a = aCoeffs;
	b = bCoeffs;
}

//*****************************************************************************
// Computees the integral over a step
//
// @param A		ODE matrix L that is used in the function
// @param yn	Current species vector solution
// @param dt	Time step to integrate over
//*****************************************************************************
VectorD BDFIntegrator::integrate(const SparseMatrixD& A, const VectorD& yn, 
	double dt){
	VectorD ynNext = 0*yn;

	// Computes the first step
	if (not initFirstStep){ computeFirstStep(A, yn, dt);};

	// Solves the BDF equation
	ynNext = solve(A, dt, previousSteps, order);

	// Shuffle the next solution into the history and removes the last
	shuffle(ynNext, previousSteps);

	return ynNext;
}

//*****************************************************************************
// Solves a BDF of a given order
//
// @param A			ODE matrix L that is used in the function
// @param dt		Time step to integrate over
// @param hist		The previous soluionts
//*****************************************************************************
VectorD BDFIntegrator::solve(const SparseMatrixD& A, double dt, MatrixD hist, 
	int BDForder){
	VectorD ynNext, rhs;
	SparseMatrixD ident(A.rows(), A.cols()), temp;
	ident.setIdentity();

	// Builds the rhs of the linear system
	rhs = buildRHS(hist, BDForder);

	// Temp vector to be inverted
	temp = ident - b(order-1)*dt*A;

	// analyze sparsisty pattern
	solver.analyzePattern(temp);

	// Compute LU decomp
	solver.factorize(temp);

	// solve
	ynNext = solver.solve(rhs);

	return ynNext;
}
//*****************************************************************************
// Computes the first step of a BDF integrator
//
// @param A			ODE matrix L that is used in the function
// @param y0	The initial condition
// @param dt		Time step to integrate over
//*****************************************************************************
void BDFIntegrator::computeFirstStep(const SparseMatrixD A, VectorD y0, 
	double dt){
	int N = 10, nTemp, tempOrder;
	double m = 1./(double)N, dtTemp;
	MatrixD tempHist; 
	VectorD sol = 0.*y0, yn = 0.*y0;

	// Sets that the first step is already calculated
	initFirstStep = true;

	// sets size of the previous steps matrix
	// needs to be the same number of rows as the initial condition and 
	// need to save order previous solutions. 
	previousSteps = MatrixD::Zero(y0.rows(), order);

	// Adds itinial condition to the history
	shuffle(y0, previousSteps);

	// This loops over lower order BDFs 
	for (int i=1; i<order; i++){
		// creats a temperary BDF order used to build the first time step
		tempOrder = i;
		// Number of steps of the temp order needed to run
		nTemp = N-i+1;
		// Time step of the lower order BDF solver
		dtTemp = dt*std::pow(m, order-i);
		// Loops over the number of times this solver needs to be called.
		// Each one of these iterations will build one of the previousSteps
		// columns. 
		std::cout << "lower bdf order " << i << std::endl;
		for (int j=0; j<=i; j++){
			// temperary history array to hold previous solutions
			tempHist = MatrixD::Zero(y0.rows(), tempOrder);
			std::cout << "index of presolve sol " << j << std::endl;
			// Loops over the number of steps for this order solver to take
			for (int step=0; step<nTemp; step++){
				// Solves over this time step
				std::cout << "index of substep " << step << std::endl;
				sol = solve(A, dtTemp, tempHist, tempOrder);
				// Adds the soluiton to temp history 
				shuffle(sol, tempHist);
			}
			nTemp = 10;
			// Adds the solution to the solvers history
			shuffle(sol, previousSteps);
		}
	}
}

//*****************************************************************************
// Addes a vector to the first colume of the history matrix and shuffels the 
// rest down a place.
//
// @param sol	The solution to be added
// @param hist	The history matrix
//*****************************************************************************
void BDFIntegrator::shuffle(const VectorD& sol, MatrixD& hist){
	int tempOrder = hist.cols();

	// Loops backward through the history matrix
	for (int i=tempOrder-1; i>0; i--){
		hist.col(i) = hist.col(i-1);
	}
	// Adds the solution vector to the first column
	hist.col(0) = sol;
}

//*****************************************************************************
// Builds the summation that is used as the right hand side of a linear solve
//
// @param hist		Matrix of previous solutions
// @param order	Order of the BDF solver
//*****************************************************************************
VectorD BDFIntegrator::buildRHS(MatrixD hist, int order){
	VectorD rhs = VectorD::Zero(hist.rows());
	VectorD aCoeff = a.row(order-1).head(order);
	// Loops over previous solutions to build the new right hand side
	for (int i=0; i < order; i++){
		rhs += aCoeff(i)*hist.col(i);
		std::cout << i << " " << aCoeff(i) << " " << hist.col(i) << std::endl;
	}

	return rhs;
}

//*****************************************************************************
// Cleans the integrator
//*****************************************************************************
void BDFIntegrator::clean(){
	initFirstStep = false;
}
//*****************************************************************************
// Method for integrator factor class
//
// @param method	explicit or implicit
// @param name		name of the solver
//*****************************************************************************
ODEintegrator *integratorFactory::getIntegrator(std::string method,
	std::string name){
	ODEintegrator *solver = nullptr;
	assert(method == "explicit" or method == "implicit");

	if (method == "explicit"){
		solver = getExplicitRKIntegrator(name);
	}
	else if (method == "implicit"){
		solver = getBDFIntegrator(name);
	}
	return solver;
}
//*****************************************************************************
// Method for integrator factor class
//
// @param name		name of the solver
//*****************************************************************************
ODEintegrator *integratorFactory::getExplicitRKIntegrator(std::string name){
	ODEintegrator *solver = nullptr;

	if (name == "forward euler"){
		// first order method
		ArrayD a(1,1); 
		ArrayD b(1,1);
		// Sets the coefficients
		a << 0.;
		b << 1.;
		solver = new explicitRKIntegrator(name, a, b);
		return solver;
	}
	else if (name == "explicit midpoint"){
		// second order method
		ArrayD a(2,2);
		ArrayD b(1,2);
		// Sets the coefficients
		a << 0.0, 0.0,
			  0.5, 0.0;
		b << 0, 1.;
		solver = new explicitRKIntegrator(name, a, b);
		return solver;
	}
	else if (name == "heun second-order"){
		// second order method
		ArrayD a(2,2);
		ArrayD b(1,2);
		// Sets the coefficients
		a << 0.0, 0.0,
			  1.0, 0.0;
		b << 0.5, 0.5;
		solver = new explicitRKIntegrator(name, a, b);
		return solver;
	}
	else if (name == "ralston second-order"){
		// second order method
		ArrayD a(2,2);
		ArrayD b(1,2);
		// Sets the coefficients
		a << 0.0, 0.0,
			  2./3., 0.0;
		b << 1./4., 3./4.;
		solver = new explicitRKIntegrator(name, a, b);
		return solver;
	}
	else if (name == "kutta third-order"){
		// third order method
		ArrayD a(3,3);
		ArrayD b(1,3);
		// Sets the coefficients
		a << 0.0, 0.0, 0.0,
			  0.5, 0.0, 0.0,
			  -1., 2.0, 0.0;
		b << 1./6., 2./3., 1./6.;
		solver = new explicitRKIntegrator(name, a, b);
		return solver;
	}
	else if (name == "heun third-order"){
		// third order method
		ArrayD a(3,3);
		ArrayD b(1,3);
		// Sets the coefficients
		a << 0.0, 0.0, 0.0,
			  1./3., 0.0, 0.0,
			  0.0, 2./3., 0.0;
		b << 1./4., 0.0, 3./4.;
		solver = new explicitRKIntegrator(name, a, b);
		return solver;
	}
	else if (name == "ralston third-order"){
		// third order method
		ArrayD a(3,3);
		ArrayD b(1,3);
		// Sets the coefficients
		a << 0.0, 0.0, 0.0,
			  0.5, 0.0, 0.0,
			  0.0, 3./4., 0.0;
		b << 2./9., 1./3., 4./9.;
		solver = new explicitRKIntegrator(name, a, b);
		return solver;
	}
	else if (name == "SSPRK3"){
		// third order method
		ArrayD a(3,3);
		ArrayD b(1,3);
		// Sets the coefficients
		a << 0.0, 0.0, 0.0,
			  1.0, 0.0, 0.0,
			  1./4., 1./4., 0.0;
		b << 1./6., 1./6., 2./3.;
		solver = new explicitRKIntegrator(name, a, b);
		return solver;
	}
	else if (name == "classic fourth-order"){
		// fourth order method
		ArrayD a(4,4);
		ArrayD b(1,4);
		// Sets the coefficients
		a << 0.0, 0.0, 0.0, 0.0,
			  1./2., 0.0, 0.0, 0.0,
			  0.0, 1./2., 0.0, 0.0,
			  0.0, 0.0, 1.0, 0.0,
		b << 1./6., 1./3., 1./3., 1./6.;
		solver = new explicitRKIntegrator(name, a, b);
		return solver;
	}
	else {
		std::string errorMessage = 
			" You have selected an explicit solver that is not\n"
			" in libowski. Avaliable solvers are:\n\n"
			" forward euler\n"
			" explicit midpoint\n"
			" heun second-order\n"
			" ralston second-order\n"
			" kutta third-order\n"
			" heun third-order\n"
			" ralston third-order\n"
			" SSPRK3\n"
			" classic fourth-order";
		libowskiException::runtimeError(errorMessage);
		return solver;
	}
}

//*****************************************************************************
// Method for integrator factor class
//
// @param name		name of the solver
//*****************************************************************************
ODEintegrator *integratorFactory::getBDFIntegrator(std::string name){
	ODEintegrator *solver = nullptr;

	if (name == "BDF1"){
		solver = new BDFIntegrator(name, 1);
		return solver;
	}
	else if (name == "BDF2"){
		solver = new BDFIntegrator(name, 2);
		return solver;
	}
	else if (name == "BDF3"){
		solver = new BDFIntegrator(name, 3);
		return solver;
	}
	else if (name == "BDF4"){
		solver = new BDFIntegrator(name, 4);
		return solver;
	}
	else if (name == "BDF5"){
		solver = new BDFIntegrator(name, 5);
		return solver;
	}
	else if (name == "BDF6"){
		solver = new BDFIntegrator(name, 6);
		return solver;
	}
	else {
		std::string errorMessage = 
			" You have selected an implicit solver\n"
			" that is not in libowski. Avaliable solvers are\n\n"
			" BDF1\n"
			" BDF2\n"
			" BDF3\n"
			" BDF4\n"
			" BDF5\n"
			" BDF6";
		libowskiException::runtimeError(errorMessage);
		return solver;
	}
}
