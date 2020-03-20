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
	assert(order > 0);

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
// Constructor for BDF integrator
//
// @param solverName		Name of the solver
// @param aCoeffs			a coefficients for BDF
// @param bCoeff			b coefficient for BDF
//*****************************************************************************
BDFIntegrator::BDFIntegrator(std::string solverName, ArrayD aCoeffs, double b):
	ODEintegrator(solverName){
	order = aCoeffs.cols();
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
	assert(order > 0);

	return ynNext;
}

//*****************************************************************************
// Computes the first step of a BDF integrator
//
// @param y0	The initial condition
//*****************************************************************************
void BDFIntegrator::computeFirstStep(VectorD y0){
	initFirstStep = true;	
}
//*****************************************************************************
// Builds the summation that is used as the right hand side of a linear solve
//
//*****************************************************************************
VectorD BDFIntegrator::buildRHS(){
	VectorD rhs;

	return rhs;
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
			" You have selected a solver that is not\n"
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
