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
#ifndef ODEINTEGRATOR_H 
#define ODEINTEGRATOR_H
#include <string>
#include <assert.h>
#include "vectorTypes.h"
#include "matrixTypes.h"
#include "exception.h"

//*****************************************************************************
// Abstract base class for ode integrator
//*****************************************************************************
class ODEintegrator{
	public:
	//**************************************************************************
	// Integrates over a single time step
	//**************************************************************************
	virtual VectorD integrate(const SparseMatrixD&, const VectorD&, double) = 0;
	//**************************************************************************
	// Constructor 
	//**************************************************************************
	ODEintegrator(std::string);
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
	// flag for if the solver is initilized
	//**************************************************************************
	bool isInit = false;
};

//*****************************************************************************
// Performs Runge-Kutta methods
//*****************************************************************************
class rungeKuttaIntegrator : public ODEintegrator{
	public:
	//**************************************************************************
	// Preforms an integration over a time step
	//**************************************************************************
	VectorD integrate(const SparseMatrixD&, const VectorD&, double) = 0;
	//**************************************************************************
	// Constructor
	//**************************************************************************
	rungeKuttaIntegrator(std::string, ArrayD, ArrayD);
	protected:
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
// Performs explicit Runge-Kutta methods
//*****************************************************************************
class explicitRKIntegrator : public rungeKuttaIntegrator{
	public:
	//**************************************************************************
	// Preforms an integration over a time step
	//**************************************************************************
	VectorD integrate(const SparseMatrixD&, const VectorD&, double);
	//**************************************************************************
	// Constructor
	//**************************************************************************
	explicitRKIntegrator(std::string, ArrayD, ArrayD);
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
	//**************************************************************************
	// Gets the ode integrator object
	//**************************************************************************
	static ODEintegrator *getIntegrator(std::string, std::string);
	private:
	//**************************************************************************
	// Gets an explicit runge integrator object
	//**************************************************************************
	static ODEintegrator *getExplicitRKIntegrator(std::string);
};
#endif
