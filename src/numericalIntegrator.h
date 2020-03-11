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
//
//	The "c" part of the Butcher tableau is used for the time dependencs on f
//	so we don't need it.
//*****************************************************************************
#ifndef NUMERICALINTEGRATOR_H 
#define NUMERICALINTEGRATOR_H
#include <string>
#include "vectorTypes.h"
#include "matrixTypes.h"
#include "exception.h"

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
	integrator(std::string, ArrayD, ArrayD);
	//**************************************************************************
	// Solver name
	//**************************************************************************
	std::string name = "None";
	protected:
	//**************************************************************************
	// The order of the method
	//**************************************************************************
	int order = -1;
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
class explicitIntegrator : public integrator{
	public:
	//**************************************************************************
	// Preforms an integration over a time step
	//**************************************************************************
	VectorD integrate(const SparseMatrixD&, const VectorD&, double);
	//**************************************************************************
	// Constructor
	//**************************************************************************
	explicitIntegrator(std::string, ArrayD, ArrayD);
	private:
	//**************************************************************************
	// The general K_i function that calculates K_i of any order 
	//**************************************************************************
	VectorD kn(int, const SparseMatrixD&, const VectorD&, double); 
};

//*****************************************************************************
// Numerical integrator factory class
//*****************************************************************************
class integratorFactory{
	public:
	static integrator *getIntegrator(std::string);
};
#endif
