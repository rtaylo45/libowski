//*****************************************************************************
// Author: Zack Taylor
// 
// Core solver for matrix expotentials. Computes y(t) = exp(A*t)v
// the solution to:
// 	y' + Ay = 0
// 
// using the CRAM Method. Returns the solution to the system.
// The code was copied from pyne CRAM solver
//*****************************************************************************
#ifndef CRAM_H
#define CRAM_H
#include <Eigen/Sparse>
#include <complex>
#include <iostream>
#include <assert.h>
#include "mpiProcess.h"
#include "linearAlgebra.h"
#include "vectorTypes.h"
#include "matrixTypes.h"


//*****************************************************************************
// The solver type; an object used to solve the matrix expotenial
//*****************************************************************************
class SolverType {

	// Private attributes
	private:
	// Poles of the radional function r
	MatrixCLD theta;
	// Residues of these poles
	MatrixCLD alpha;
	// Limit of r at infinity
	long double alpha_0 = 0.0L;

	//*************************************************************************
	// Initialization of solver
	//*************************************************************************
	public:
	SolverType();
	//*************************************************************************
	// Solver function
	//*************************************************************************
	VectorD solve(const SparseMatrixD&, const VectorD&, double);

	//*************************************************************************
	// Sets the solver type for CRAM
	//*************************************************************************
	void setSolveType(std::string);

	private:
	//*************************************************************************
	// Calculates the quadrature points for the parabolic contour
	//*************************************************************************
	MatrixCLD parabolicContourCoeffs(int);

	//*************************************************************************
	// Calculates the quadrature points for the hyperbolic contour
	//*************************************************************************
	MatrixCLD hyperbolicContourCoeffs(int);

	//*************************************************************************
	// Solver pointer to the method of matrix exponential solve type. 
	// Default is the base solve without scaling and squaring.
	//*************************************************************************
	VectorD (SolverType::*solverPtr)(const SparseMatrixD&, const VectorD&, 
		double) = &SolverType::solveBase;

	//*************************************************************************
	// Solves CRAM with no matrix scaling
	//*************************************************************************
	VectorD solveBase(const SparseMatrixD&, const VectorD&, double);

	//*************************************************************************
	// Solves CRAM with matrix scaling
	//*************************************************************************
	VectorD solveScale(const SparseMatrixD&, const VectorD&, double);
};
#endif
