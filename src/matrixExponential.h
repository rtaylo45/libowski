//*****************************************************************************
// Author: Zack Taylor
// 
// Matrix exponential class that defines the various methods to solve the 
// matrix exponential. The solver can run two mthods, apply and compute.
// Apply, directly computes the action of the matrix exponential on a vector.
// Compute, directly computes the matrix expontial its self. Use apply when
// possible over compute because compute requires the inverse of a matrix
//*****************************************************************************
#ifndef MATRIXEXPONENTIAL_H
#define MATRIXEXPONENTIAL_H
#include <Eigen/Sparse>
#include <complex>
#include <string>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include "mpiProcess.h"
#include "linearAlgebra.h"
#include "vectorTypes.h"
#include "matrixTypes.h"
#include "exception.h"
#include "utilBase.h"

//*****************************************************************************
// Abstract base matrix exponential class
//*****************************************************************************
class matrixExponential{
	public:
	//**************************************************************************
	// Computes the matrix exponential action on a vector. exp(A*t)v
	//**************************************************************************
	virtual VectorD apply(const SparseMatrixD&, const VectorD&, double) = 0;
	//**************************************************************************
	// Computes the matrix exponential. exp(A*t)
	//**************************************************************************
	virtual SparseMatrixD compute(const SparseMatrixD&, double) = 0;
	//**************************************************************************
	// Consturcture
	//**************************************************************************
	matrixExponential(bool, int);	
	//**************************************************************************
	// Solver name
	//**************************************************************************
	std::string name = "None";
	protected:
	//**************************************************************************
	// Logic to set if the krylov subspace method should be applied to the 
	// solver
	//**************************************************************************
	bool useKrylovSubspace = false;
	//**************************************************************************
	// Krylov subspace dimension
	//**************************************************************************
	int krylovSubspaceDim = 10;

};

//*****************************************************************************
// Abstract base class for matrix exponential methods based on the pade
// approximation
//*****************************************************************************
class pade : public matrixExponential{
	protected:
	//**************************************************************************
	// Pade approximetion of order (3,3) 
	//**************************************************************************
	void pade3(const SparseMatrixD&, const SparseMatrixD&, SparseMatrixD&, 
		SparseMatrixD&);
	//**************************************************************************
	// Pade approximetion of order (5,5) 
	//**************************************************************************
	void pade5(const SparseMatrixD&, const SparseMatrixD&, const SparseMatrixD&,
		SparseMatrixD&, SparseMatrixD&);
	//**************************************************************************
	// Pade approximetion of order (7,7) 
	//**************************************************************************
	void pade7(const SparseMatrixD&, const SparseMatrixD&, const SparseMatrixD&,
		const SparseMatrixD&, SparseMatrixD&, SparseMatrixD&);
	//**************************************************************************
	// Pade approximetion of order (9,9) 
	//**************************************************************************
	void pade9(const SparseMatrixD&, const SparseMatrixD&, const SparseMatrixD&,
		const SparseMatrixD&, const SparseMatrixD&, SparseMatrixD&, SparseMatrixD&);
	//**************************************************************************
	// Pade approximetion of order (13,13) 
	//**************************************************************************
	void pade13(const SparseMatrixD&, const SparseMatrixD&, const SparseMatrixD&,
		const SparseMatrixD&, SparseMatrixD&, SparseMatrixD&);
	//**************************************************************************
	// Runs the algorithm
	//**************************************************************************
	virtual void run(const SparseMatrixD&, SparseMatrixD&, SparseMatrixD&, int&)=0;

	public:
	//**************************************************************************
	// Constructor
	//**************************************************************************
	pade(bool, int);
	//**************************************************************************
	// Computes the matrix exponential action on a vector. exp(A*t)v
	//**************************************************************************
	virtual VectorD apply(const SparseMatrixD&, const VectorD&, double);
	//**************************************************************************
	// Computes the matrix exponential. exp(A*t)
	//**************************************************************************
	virtual SparseMatrixD compute(const SparseMatrixD&, double);
};

//*****************************************************************************
// Pade class for method 1. Based on the algorithm
//
// The scaling and squaring method for the matrix exponential
// revisited.
// SIAM Journal on  Matrix Analysis  and  Applications, 26(4):1179–1193, 2005
//*****************************************************************************
class method1 : public pade{
	public:
	//**************************************************************************
	// Constructor
	//**************************************************************************
	method1(bool, int);
	private:
	//**************************************************************************
	// Runs the algorithm
	//**************************************************************************
	virtual void run(const SparseMatrixD&, SparseMatrixD&, SparseMatrixD&, int&);	
};

//*****************************************************************************
// Pade class for method 2. Based on the algorithm
//
// A new scaling and squaring algorithm for the matrix exponential.
// SIAM  Journal  on  Matrix  Analysis  and  Applications, 31, 01 2009
//*****************************************************************************
class method2 : public pade{
	public:
	//**************************************************************************
	// Constructor
	//**************************************************************************
	method2(bool, int);
	
	private:
	//**************************************************************************
	// Runs the algorithm
	//**************************************************************************
	virtual void run(const SparseMatrixD&, SparseMatrixD&, SparseMatrixD&, int&);
	
	//**************************************************************************
	// Normest, some fucntion they define in updated method 
	//**************************************************************************
	double normest(const SparseMatrixD&, const SparseMatrixD&);
	
	//**************************************************************************
	// Normest, some fucntion they define in updated method 
	//**************************************************************************
	double normest(const SparseMatrixD&, const int);

	//**************************************************************************
	// ell, some fucntion they define in updated method 
	//**************************************************************************
	int ell(const SparseMatrixD&, const int);
};

//*****************************************************************************
// Abstract base class for matrix exponential methods based on Cauchys
// integral formula
//*****************************************************************************
class cauchy : public matrixExponential{
	protected:
	// Poles of the rational function r
	MatrixCLD theta;
	// Residues of these poles
	MatrixCLD alpha;
	// Limit of r at infinity
	long double alpha_0 = 0.0L;

	public:
	//**************************************************************************
	// Constructer
	//**************************************************************************
	cauchy(bool, int);
	//**************************************************************************
	// Computes the matrix exponential action on a vector. exp(A*t)v
	//**************************************************************************
	virtual VectorD apply(const SparseMatrixD&, const VectorD&, double);
	//**************************************************************************
	// Computes the matrix exponential. exp(A*t)
	//**************************************************************************
	virtual SparseMatrixD compute(const SparseMatrixD&, double);
};

//*****************************************************************************
// Implimentation of the CRAM algorithm for solving the matrix exponential.
// Currently the only order avaliable is 16.
//*****************************************************************************
class CRAM : public cauchy{
	public:
	//**************************************************************************
	// Constructor for the CRAM class
	//**************************************************************************
	CRAM(bool krylovBool, int krylovDim);
};

//*****************************************************************************
// Implimentation Cauchys integral formula with a parabolic contour around
// the negative real axis. Coefficients for any order can be computed but
// 32 is the default.
//*****************************************************************************
class parabolic : public cauchy{
	public:
	//**************************************************************************
	// Constructor for the parabolic class
	//**************************************************************************
	parabolic(bool krylovBool, int krylovDim);

	private:
	// Order of the approximation
	int order = 32;
	//*************************************************************************
	// Calculates the quadrature points for the parabolic contour
	//*************************************************************************
	MatrixCLD parabolicContourCoeffs(int);
};

//*****************************************************************************
// Implimentation Cauchys integral formula with a hyperbolic contour around
// the negative real axis. Coefficients for any order can be computed but
// 32 is the default.
//*****************************************************************************
class hyperbolic : public cauchy{
	public:
	//**************************************************************************
	// Constructor for the hyperbolic class
	//**************************************************************************
	hyperbolic(bool krylovBool, int krylovDim);

	private:
	// Order of the approximation
	int order = 32;
	//**************************************************************************
	// Calculates the quadrature points for the hyperbolic contour
	//**************************************************************************
	MatrixCLD hyperbolicContourCoeffs(int);
};
//*****************************************************************************
// Matrix exponential factory class
//*****************************************************************************
class matrixExponentialFactory{
	public:
	static matrixExponential *getExpSolver(std::string, bool = false, int = 10);
};
#endif 
