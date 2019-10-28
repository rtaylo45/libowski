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
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <complex>
#include <iostream>
#include "mpiProcess.h"
#include "linearAlgebra.h"


//*****************************************************************************
// The solver type; an object used to solve the matrix expotenial
//*****************************************************************************
class SolverType {

	// Private attributes
	private:
	// Poles of the radional function r
	Eigen::Matrix<std::complex<double>,8,1> theta;
	// Residues of these poles
	Eigen::Matrix<std::complex<double>,8,1> alpha;
	// Limit of r at infinity
	double alpha_0 = 2.1248537104952237488e-16;

	//*************************************************************************
	// Initialization of solver
	//*************************************************************************
	public:
	SolverType();
	//*************************************************************************
	// Solver function
	//*************************************************************************
	Eigen::MatrixXd solve(Eigen::SparseMatrix<double>, Eigen::VectorXd, double);
	private:
	//*************************************************************************
	// Builds a sparse identity matrix
	//*************************************************************************
	Eigen::SparseMatrix<std::complex<double>> buildSparseIdentity(int n);

	//*************************************************************************
	// Solves CRAM with no matrix scaling
	//*************************************************************************
	Eigen::MatrixXd solveBase(Eigen::SparseMatrix<double>, Eigen::VectorXd, double);

	//*************************************************************************
	// Solves CRAM with matrix scaling
	//*************************************************************************
	Eigen::MatrixXd solveScale(Eigen::SparseMatrix<double>, Eigen::VectorXd, double);
};
#endif
