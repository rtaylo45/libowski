//*****************************************************************************
// Author: Zack Taylor
// 
// Provies a class of numerical ode integrators 
// Classical numerical integrator class that implements ode integrator methods
// These methods solve the general form,
//
//		dy/dt = f(t,y).
//
//	In libowski this takes the form,
//
//		dy/dt = Ly.
//
//	If f(t,y) = Ly then the function is not dependent on t and 
//	f(t,y) = f(y). 
//*****************************************************************************
#ifndef ODEINTEGRATOR_H 
#define ODEINTEGRATOR_H
#include <string>
#include <assert.h>
#include <math.h>
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
	// Cleans the solver
	//**************************************************************************
	virtual void clean()=0;
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
	//**************************************************************************
	// Cleans the solver
	//**************************************************************************
	void clean()=0;
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
	//**************************************************************************
	// Cleans the solver
	//**************************************************************************
	void clean();
	private:
	//**************************************************************************
	// The general K_i function that calculates K_i of any order 
	//**************************************************************************
	VectorD kn(int, const SparseMatrixD&, const VectorD&, double); 
};

//*****************************************************************************
// Preforms backward differentiation method. This is a class of implicit 
// solvers
//
// The general formula ofr a BDF can be written as
//
//		Sum a_k*y_(n+k) = h*b*L*y_(n+s)	from k = 0 to k = s
//
// where s is the order of the method. Solving for y_(n+s)
//
//		y_(n+s) = (I - b*h*L)^(-1) * (Sum a_k*y_(n+k))  from k = 0 to k = s - 1
//*****************************************************************************
class BDFIntegrator : public ODEintegrator{
	public:
	//**************************************************************************
	// Preforms an integration over a time step
	//**************************************************************************
	VectorD integrate(const SparseMatrixD&, const VectorD&, double);
	//**************************************************************************
	// Constructor
	//**************************************************************************
	BDFIntegrator(std::string, int);
	//**************************************************************************
	// Cleans the solver
	//**************************************************************************
	void clean();
	private:
	//**************************************************************************
	// Solves the BDF equation
	//**************************************************************************
	VectorD solve(const SparseMatrixD&, double, MatrixD, int);
	//**************************************************************************
	// Addes a solution vector to the first column of the pervious soluitons
	// and moves the rest back
	//**************************************************************************
	void shuffle(const VectorD&, MatrixD&);
	//**************************************************************************
	// computes the first step
	//**************************************************************************
	void computeFirstStep(SparseMatrixD, VectorD, double);
	//**************************************************************************
	// Builds the rhs vector for the linear implicit solve
	//**************************************************************************
	VectorD buildRHS(MatrixD, int);
	//**************************************************************************
	// logical stating if the first step has been computed or not
	//**************************************************************************
	bool initFirstStep = false;
	//**************************************************************************
	// Matrix containing previous solutions
	//**************************************************************************
	MatrixD previousSteps;
	//**************************************************************************
	// a coefficients for BDF. Each order is house in a row of the array, with
	// the oppisite sign
	//**************************************************************************
	MatrixD a;

	//**************************************************************************
	// b coefficients for BDF.
	//**************************************************************************
	VectorD b;
	//**************************************************************************
	// Sparse LU Linear Solver
	//**************************************************************************
	Eigen::SparseLU<SparseMatrixD, COLAMDOrdering<int> > solver;
	
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
	// Gets an explicit runge kutta integrator object
	//**************************************************************************
	static ODEintegrator *getExplicitRKIntegrator(std::string);
	//**************************************************************************
	// Gets a BDF integrator object
	//**************************************************************************
	static ODEintegrator *getBDFIntegrator(std::string);
};
#endif
