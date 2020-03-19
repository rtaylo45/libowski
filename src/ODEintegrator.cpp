//*****************************************************************************
// Author: Zack Taylor
// 
// Classical numerical integrator class that implements Runge-Kutta
// methods. These methods solve the general form,
//*****************************************************************************
#include "ODEintegrator.h"

//*****************************************************************************
// Default constructor for ode integrator class
//*****************************************************************************
ODEintegrator::ODEintegrator(std::string solverName){
	name = solverName;
	isInit = true;
}

//*****************************************************************************
// Constructor for explicit integrator class
//*****************************************************************************
rungeKuttaIntegrator::rungeKuttaIntegrator(std::string solverName, ArrayD aCoeffs,
	ArrayD bCoeffs):ODEintegrator(solverName){
	order = aCoeffs.cols();
	a = aCoeffs;
	b = bCoeffs;
}
//*****************************************************************************
// Methods for the explicit integrator
//
// @param A		ODE matrix L that is used in the function
// @param yn	Current species vector solution
// @param dt	Time step to integrate over
//*****************************************************************************
VectorD rungeKuttaIntegrator::integrate(const SparseMatrixD& A, const VectorD& yn,
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
VectorD rungeKuttaIntegrator::kn(int order, const SparseMatrixD& A, const
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
// Method for integrator factor class
//*****************************************************************************
ODEintegrator *integratorFactory::getIntegrator(std::string type){
	ODEintegrator *solver = nullptr;

	if (type == "forward euler"){
		// first order method
		ArrayD a(1,1); 
		ArrayD b(1,1);
		// Sets the coefficients
		a << 0.;
		b << 1.;
		solver = new rungeKuttaIntegrator(type, a, b);
		return solver;
	}
	else if (type == "explicit midpoint"){
		// second order method
		ArrayD a(2,2);
		ArrayD b(1,2);
		// Sets the coefficients
		a << 0.0, 0.0,
			  0.5, 0.0;
		b << 0, 1.;
		solver = new rungeKuttaIntegrator(type, a, b);
		return solver;
	}
	else if (type == "heun second-order"){
		// second order method
		ArrayD a(2,2);
		ArrayD b(1,2);
		// Sets the coefficients
		a << 0.0, 0.0,
			  1.0, 0.0;
		b << 0.5, 0.5;
		solver = new rungeKuttaIntegrator(type, a, b);
		return solver;
	}
	else if (type == "ralston second-order"){
		// second order method
		ArrayD a(2,2);
		ArrayD b(1,2);
		// Sets the coefficients
		a << 0.0, 0.0,
			  2./3., 0.0;
		b << 1./4., 3./4.;
		solver = new rungeKuttaIntegrator(type, a, b);
		return solver;
	}
	else if (type == "kutta third-order"){
		// third order method
		ArrayD a(3,3);
		ArrayD b(1,3);
		// Sets the coefficients
		a << 0.0, 0.0, 0.0,
			  0.5, 0.0, 0.0,
			  -1., 2.0, 0.0;
		b << 1./6., 2./3., 1./6.;
		solver = new rungeKuttaIntegrator(type, a, b);
		return solver;
	}
	else if (type == "heun third-order"){
		// third order method
		ArrayD a(3,3);
		ArrayD b(1,3);
		// Sets the coefficients
		a << 0.0, 0.0, 0.0,
			  1./3., 0.0, 0.0,
			  0.0, 2./3., 0.0;
		b << 1./4., 0.0, 3./4.;
		solver = new rungeKuttaIntegrator(type, a, b);
		return solver;
	}
	else if (type == "ralston third-order"){
		// third order method
		ArrayD a(3,3);
		ArrayD b(1,3);
		// Sets the coefficients
		a << 0.0, 0.0, 0.0,
			  0.5, 0.0, 0.0,
			  0.0, 3./4., 0.0;
		b << 2./9., 1./3., 4./9.;
		solver = new rungeKuttaIntegrator(type, a, b);
		return solver;
	}
	else if (type == "SSPRK3"){
		// third order method
		ArrayD a(3,3);
		ArrayD b(1,3);
		// Sets the coefficients
		a << 0.0, 0.0, 0.0,
			  1.0, 0.0, 0.0,
			  1./4., 1./4., 0.0;
		b << 1./6., 1./6., 2./3.;
		solver = new rungeKuttaIntegrator(type, a, b);
		return solver;
	}
	else if (type == "classic fourth-order"){
		// fourth order method
		ArrayD a(4,4);
		ArrayD b(1,4);
		// Sets the coefficients
		a << 0.0, 0.0, 0.0, 0.0,
			  1./2., 0.0, 0.0, 0.0,
			  0.0, 1./2., 0.0, 0.0,
			  0.0, 0.0, 1.0, 0.0,
		b << 1./6., 1./3., 1./3., 1./6.;
		solver = new rungeKuttaIntegrator(type, a, b);
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
