//*****************************************************************************
// Author: Zack Taylor
// 
// Classical numerical integrator class that implements explicit Runge-Kutta
// methods. These methods solve the general form,
//*****************************************************************************
#include "numericalIntegrator.h"

//*****************************************************************************
// Default constructor for numerical integrator class
//*****************************************************************************
integrator::integrator(std::string solverName, ArrayD aCoeffs, ArrayD bCoeffs){
	name = solverName;
	a = aCoeffs;
	b = bCoeffs;
}

//*****************************************************************************
// Constructor for explicit integrator class
//*****************************************************************************
explicitIntegrator::explicitIntegrator(std::string solverName, ArrayD aCoeffs,
	ArrayD bCoeffs):integrator(solverName, aCoeffs, bCoeffs){};
//*****************************************************************************
// Methods for the explicit integrator
//
// @param A		ODE matrix L that is used in the function
// @param yn	Current species vector solution
// @param dt	Time step to integrate over
//*****************************************************************************
VectorD explicitIntegrator::integrate(const SparseMatrixD& A, const VectorD& yn,
	double dt){
	VectorD ynPlus1 = 0*yn, ySum = 0*yn;

	// Loops over the order	
	for (int i=1; i <= order; i++){
		ySum = b[i-1]*k(i, A, yn, dt);	
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
VectorD explicitIntegrator::kn(int order, const SparseMatrixD& A, const
	VectorD& yn, double dt){
	VectorD sum = 0*yn, tempy = 0*yn;

	// Sums over kn's to build the current kn	
	for (int j=1; j < order; j++){
		sum += a[order-1,j-1]*kn(j, A, yn, dt);
	}
	tempy = yn + dt*sum;
	return A*tempy;
}
