//*****************************************************************************
// Author: Zack Taylor
// 
// Classical numerical integrator class that implements explicit Runge-Kutta
// methods. These methods solve the general form,
//
//		dy/dt = f(t,y).
//
//	In libowski this takes the form,
//
//		dy/dt = Ly.
//
//	If f(t,y) = Ly then the function is not dependent on t and 
//	f(t,y) = f(y). 
//
//	Each Runge-Kutta method follows the same general formula,
//
//		y_(n+1) = y_n + h*Sum (b_i * k_i)	from i-1 to s
//
//		k_i = f[y_n + h*Sum (a_(i,j)*k_j)]	from j=1 to i-1
//
//	The coefficients take the form of the Butcher tableau
//
//		c_1 | a_(1,1)	a_(1,2) ...  a_(1,s)
//		c_2 | a_(2,1)	a_(2,2) ...	 a_(2,s)
//		 .  |    .        .            .
//		 .  |    .        .            .
//		 .  |    .        .            .
//		c_s | a_{s,1)  a_(2,s) ...  a_(s,s)
//*****************************************************************************
#ifndef NUMERICALINTEGRATOR_H 
#define NUMERICALINTEGRATOR_H
#include "vectorTypes.h"
#include "matrixTypes.h"
#include <string>

//*****************************************************************************
// Abstract base class for numerical integrator
//*****************************************************************************
class integrator{
	public:
	//**************************************************************************
	// Integrates over a single time step
	//**************************************************************************
	virtual VectorD integrate(const SparseMatrixD&, const VectorD&, double) = 0;
	//**************************************************************************
	// Constructor 
	//**************************************************************************
	numericalIntegrator(std::string);
	//**************************************************************************
	// Solver name
	//**************************************************************************
	std::sring name = "None";
	protected:
	//**************************************************************************
	// The order of the method
	//**************************************************************************
	int order;
	//**************************************************************************
	// Holds the a coefficients for k functions
	//**************************************************************************
	ArrayD a;
	//**************************************************************************
	// Holds the b coefficients for Runge-Kutta methods
	//**************************************************************************
	ArrayD b;
};

//*****************************************************************************
// Explitict integrator class. Performs explicit Runge-Kutta methods
//*****************************************************************************
class integrator : public explicitIntegrator{
	//**************************************************************************
	// Preforms an integration over a time step
	//**************************************************************************
	VectorD integrate(const SparseMatrixD&, const VectorD&, double);
	private:
	//**************************************************************************
	// The general K_i function that calculates K_i of any order 
	//**************************************************************************
	double kn(int, const SparseMatrixD&, const VectorD&, double); 
	//**************************************************************************
	// Computes k of order 1
	//**************************************************************************
	double k1(const SparseMatrixD&, const VectorD&);
	//**************************************************************************
	// Computes k of order 2
	//**************************************************************************
	double k2(const SparseMatrixD&, const VectorD&, double);
	//**************************************************************************
	// Computes k of order 3
	//**************************************************************************
	double k3(const SparseMatrixD&, const VectorD&, double);
	//**************************************************************************
	// Computes k of order 4
	//**************************************************************************
	double k4(const SparseMatrixD&, const VectorD&, double);
};
#endif
