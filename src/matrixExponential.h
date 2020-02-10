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
#include "mpiProcess.h"
#include "linearAlgebra.h"
#include "vectorTypes.h"
#include "matrixTypes.h"

//*****************************************************************************
// Abstract base matrix exponential class
//*****************************************************************************
class matrixExponential{
	public:
	//*****************************************************************************
	// Computes the matrix exponential action on a vector. exp(A*t)v
	//*****************************************************************************
	virtual VectorD apply(const SparseMatrixD&, const VectorD&, double) = 0;
	//*****************************************************************************
	// Computes the matrix exponential. exp(A*t)
	//*****************************************************************************
	virtual SparseMatrixD compute(const SparseMatrixD&, double) = 0;
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
	//*****************************************************************************
	// Computes the matrix exponential action on a vector. exp(A*t)v
	//*****************************************************************************
	virtual VectorD apply(const SparseMatrixD&, const VectorD&, double);
	//*****************************************************************************
	// Computes the matrix exponential. exp(A*t)
	//*****************************************************************************
	virtual SparseMatrixD compute(const SparseMatrixD&, double);
};

//*****************************************************************************
// Implimentation of the CRAM algorithm for solving the matrix exponential.
// Currently the only order avaliable is 16.
//*****************************************************************************
class CRAM : public cauchy{
	public:
	//*****************************************************************************
	// Constructor for the CRAM class
	//*****************************************************************************
	CRAM();
};

//*****************************************************************************
// Implimentation Cauchys integral formula with a parabolic contour around
// the negative real axis. Coefficients for any order can be computed but
// 32 is the default.
//*****************************************************************************
class parabolic : public cauchy{
	public:
	//*****************************************************************************
	// Constructor for the parabolic class
	//*****************************************************************************
	parabolic();

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
	//*****************************************************************************
	// Constructor for the hyperbolic class
	//*****************************************************************************
	hyperbolic();

	private:
	// Order of the approximation
	int order = 32;
	//*************************************************************************
	// Calculates the quadrature points for the hyperbolic contour
	//*************************************************************************
	MatrixCLD hyperbolicContourCoeffs(int);
};
//*****************************************************************************
// Matrix exponential factory class
//*****************************************************************************
class matrixExponentialFactory{
	public:
	static matrixExponential *getExpSolver(std::string type);
};
#endif 
