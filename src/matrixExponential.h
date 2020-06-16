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
#include <fstream>
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
	//**************************************************************************
	// Logic to se if you want to balance the matrix (can help to reduce the 
	// matrix norm... use with caution) This is needed to methods whoes 
	// accuracy deponds on the matrix norm.
	//**************************************************************************
	bool useBalance = false;
};

//*****************************************************************************
// Class for computing the action of a matrix exponential on a vector using 
// a Taylor series expansion from the following paper:
//
// Al-Mohy, Awad H. and Higham, Nicholas J. (2011) Computing the Action of 
// the Matrix Exponential, with an Application to Exponential Integrators. 
// SIAM Journal on Scientific Computing, 33 (2). pp. 488-511. ISSN 1064-8275
//*****************************************************************************
class taylor : public matrixExponential{
	private:
	// Theta
	VectorD theta = VectorD::Zero(100);
	//**************************************************************************
	// Selects the Taylor series degree for the approximation
	//**************************************************************************
	void parameters(const SparseMatrixD&, const VectorD&, MatrixD&, int = 55, 
		int = 8, bool = true, bool = false);
	//**************************************************************************
	// Internal fucntion that computes the action of a matrix exponential on 
	// a vector. exp(A*t)v
	//**************************************************************************
	VectorD expmv(const SparseMatrixD&, const double, const VectorD&, MatrixD&,
		bool = true, bool = false);


	public:
	//**************************************************************************
	// Constructor
	//**************************************************************************
	taylor(bool, int, bool);
	//**************************************************************************
	// Computes the matrix exponential action on a vector. exp(A*t)v
	//**************************************************************************
	VectorD apply(const SparseMatrixD&, const VectorD&, double);
	//**************************************************************************
	// Computes the matrix exponential. exp(A*t)
	//**************************************************************************
	SparseMatrixD compute(const SparseMatrixD&, double);
	
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
	pade(bool, int, bool);
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
	method1(bool, int, bool);
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
	method2(bool, int, bool);
	
	private:
	//**************************************************************************
	// Runs the algorithm
	//**************************************************************************
	virtual void run(const SparseMatrixD&, SparseMatrixD&, SparseMatrixD&, int&);
	
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
	cauchy(bool, int, bool);
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
	CRAM(bool, int, bool);
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
	parabolic(bool, int, bool);

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
	hyperbolic(bool, int, bool);

	private:
	// Order of the approximation
	int order = 32;
	//**************************************************************************
	// Calculates the quadrature points for the hyperbolic contour
	//**************************************************************************
	MatrixCLD hyperbolicContourCoeffs(int);
};

//*****************************************************************************
// Computes the matrix exponential using Laguerre Polynomials 
//
// "Using Generalized Laguerre Polynomials to Compute the Matrix
// Exponential in Burnup Equations" - Ding She
//*****************************************************************************
class LPAM : public matrixExponential{
	private:
	// Sparse LU solver
	Eigen::SparseLU<SparseMatrixD, COLAMDOrdering<int>> solver;
	// Matrix inverse
	SparseMatrixD matInverse;
	// Scaling factor for Laguerre Polynomials
	double tau = 20.;
	// Some parameter for Laguerre Polynomials
	double a = 20.;
	// Number of expansion polynomials
	int k = 18;
	//**************************************************************************
	// Calculates the leading Laguerre coefficient
	//**************************************************************************
	double laguerreCoefficient(double, double, double);
	//**************************************************************************
	// Calculates the polynomal for apply
	//**************************************************************************
	VectorD laguerrePolynomial(double, int, int, const VectorD&, const 
		SparseMatrixD&);
	//**************************************************************************
	// Calculates the polynomal for compute
	//**************************************************************************
	SparseMatrixD laguerrePolynomial(double, int, int, const SparseMatrixD&);

	public:
	//**************************************************************************
	// Constructor
	//**************************************************************************
	LPAM(bool, int, bool);
	//**************************************************************************
	// Computes the matrix exponential action on a vector. exp(A*t)v
	//**************************************************************************
	VectorD apply(const SparseMatrixD&, const VectorD&, double);
	//**************************************************************************
	// Computes the matrix exponential. exp(A*t)
	//**************************************************************************
	SparseMatrixD compute(const SparseMatrixD&, double);
};

//*****************************************************************************
// Matrix exponential factory class
//*****************************************************************************
class matrixExponentialFactory{
	public:
	static matrixExponential *getExpSolver(std::string, bool = false, int = 10, 
		bool = false);
};
#endif 
