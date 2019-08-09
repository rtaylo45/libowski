#include <Eigen/Core>
#include "coreSolver.h"
#include <chrono>
#include <assert.h>
#include <iostream>
#include <Eigen/Sparse>
#include <vector>
#include <random>

using namespace std::chrono;
using namespace Eigen;

bool isApprox(double val1, double val2, double tol = 1e-10){

	double diff = abs(val1 - val2);
	if (diff < tol) { return true; }
	else { return false; }
}

void testSolverTime(){
//*****************************************************************************
//	Problem statement:
//		A is a 100,000 x 100,000 size tridiagonal matrix. The matrix entries
//		are random numbers between 0 and 1. The problem is ran 20 times to 
//		ensure each solver time is below 0.8 seconds.
//
//*****************************************************************************
	typedef Eigen::Triplet<double> T;
	std::default_random_engine gen;
	std::uniform_real_distribution<double> dist(0.0,1.0);

	for (int i = 1; i < 20; i++){
		int n = 100000;

		// Get the solver object
		SolverType ExpSolver;

		// Build random dense matrix
		//MatrixXd A = MatrixXd::Random(n,n);
		//MatrixXd n0 = MatrixXd::Random(n,1);

    	SparseMatrix<double> A(n,n);
    	SparseVector<double> n0(n);
		MatrixXd sol;
		std::vector<T> tripletList;
		tripletList.reserve(3*n);

		for (int j = 0; j < n; j++){
			auto v_ij=dist(gen); //generate random number
			if(v_ij > 0.001){ tripletList.push_back(T(j,j,v_ij)); };

			if(j < n-2)
			{
				auto v_ij=dist(gen); //generate random number
				if(v_ij > 0.001){ tripletList.push_back(T(j,j+2,v_ij)); };
			}
			if(j > 1)
			{
				auto v_ij=dist(gen); //generate random number
				if(v_ij > 0.001){ tripletList.push_back(T(j,j-2,v_ij)); };
			}

		}

		A.setFromTriplets(tripletList.begin(), tripletList.end());

		// Start time
		auto start = high_resolution_clock::now();
		sol = ExpSolver.solve(A, n0, 1.0);
		auto end = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(end - start);

		// Convert to seconds
		assert(duration.count()/1.e6 < 0.8);
		//std::cout << "Time: " << duration.count()/1.e6 << std::endl;

	}
}

void tankProblem(){
//*****************************************************************************
//	Problem statement:
//		Let brine tanks 1, 2, 3 be given of volumes 20, 40, 60, It is supposed 
//		that fluid enters tank A at rate r, drains from A to B at rate r, 
//		drains from B to C at rate r, then drains from tank C at rate r. Hence 
//		the volumes of the tanks remain constant. Let r = 10. The problem is 
//		Let x1(t), x2(t), x3(t) denote the amount of salt at time t in each 
//		tank. the problem is taken from:
//		www.math.utah.edu/~gustafso/2250systems-de-1.pdf
//
//	Initial conditons:
//		x1_0 = 1000.0
//		x2_0 = 0.0
//		x3_0 = 0.0
//
//*****************************************************************************
	typedef Eigen::Triplet<double> T;
    double x1_0 = 1000.0, x2_0 = 0.0, x3_0 = 0.0;
    double t = 0.0; 
	double x1, x2, x3;
	int steps = 100;
	double totalTime = 50.0;
	double dt = totalTime/steps;
    SparseMatrix<double> A(3,3);
    SparseVector<double> b(3);
	MatrixXd sol;
	std::vector<T> tripletList;
	tripletList.reserve(5);

	// Set the A matrix and inition condition vector b
	b.insert(0,0) = x1_0;
	tripletList.push_back(T(0,0,-0.5)); tripletList.push_back(T(1,0,0.5));
	tripletList.push_back(T(1,1,-0.25)); tripletList.push_back(T(2,1,0.25));
	tripletList.push_back(T(2,2,-1./6.)); 
	A.setFromTriplets(tripletList.begin(), tripletList.end());

	// Sets the solver
    SolverType ExpSolver;

	for (int i = 0; i < steps; i++){

		t = t + dt;

		sol = ExpSolver.solve(A, b, t);
    	x1 = x1_0*exp(-t/2.);
    	x2 = -2.*x1_0*exp(-t/2.) + (x2_0 + 2.*x1_0)*exp(-t/4.);
    	x3 = (3./2.)*x1_0*exp(-t/2.) - 3.*(x3_0 + 2.*x1_0)*exp(-t/4) +
    	    (x3_0 - (3./2.)*x1_0 + 3.*(x2_0 + 2.*x1_0))*exp(-t/6.);
	
    	std::cout << x1 << " " << sol(0) << std::endl;
    	std::cout << x2 << " " << sol(1) << std::endl;
    	std::cout << x3 << " " << sol(2) << std::endl;
    	std::cout << " " << std::endl;

		assert(isApprox(sol(0), x1));
		assert(isApprox(sol(1), x2));
		assert(isApprox(sol(2), x3));
	}
}

void xenonIodineProblem(){
//*****************************************************************************
//	Problem statement:
//		dN_xe/dt = gamma_xe*Sigma_f*flux - sigma_a*flux*N_xe + lamba_I*N_I 
//			- lambda_xe*N_xe
//
//		dN_I/dt = gamma_I*Sigma_f*flux - lambda_I*N_I
//
//	Initial conditons:
//		N_xe_0 = 0.0
//		N_I_0 = 0.0
//
//	To add the constant source terms we need to add a dummy species to hold the
//	coefficients. 
//		
//		dN_d/dt = 0.0
//		d_0 = 1.0
//	
//*****************************************************************************

	typedef Eigen::Triplet<double> T;
    double N_xe_0 = 0.0, N_I_0 = 0.0, N_d_0 = 1.0;
    double t = 0.0; 
	double N_xe, N_I, N_d;
	double a, b, c;
	int steps = 1;
	double totalTime = 100.0;
	double dt = totalTime/steps;
	double lambda_I = 2.11E-5;
	double lambda_xe = 2.9306E-5;
	double sigma_a = 2.002E-22;
	double Sigma_f = 9.7532E-1;
	double flux = 2.5E16;
	double gamma_xe = 0.002468;
	double gamma_I = 0.063033;
    SparseMatrix<double> A(3,3);
    SparseVector<double> N0(3);
	MatrixXd sol;
	std::vector<T> tripletList;
	tripletList.reserve(5);

	// Sets the solver
    SolverType ExpSolver;

    //A = MatrixXd::Zero(3,3);
    //N0 = MatrixXd::Zero(3,1);

	//A(0,0) = lambda_I; A(0,1) = lambda_I - sigma_a*flux;
	//A(0,2) = gamma_xe*Sigma_f*flux; A(1,1) = -lambda_I;
	//A(1,2) = gamma_I*Sigma_f*flux; N0(2,0) = N_d_0;

	N0.insert(2,0) = N_d_0;
	tripletList.push_back(T(0,0,lambda_I)); 
	tripletList.push_back(T(0,1,lambda_I - sigma_a*flux));
	tripletList.push_back(T(0,2,gamma_xe*Sigma_f*flux)); 
	tripletList.push_back(T(1,1,-lambda_I));
	tripletList.push_back(T(1,2,gamma_I*Sigma_f*flux)); 
	A.setFromTriplets(tripletList.begin(), tripletList.end());

	for (int i = 0; i < steps; i++){

		t = t + dt;

		sol = ExpSolver.solve(A, N0, t);

		a = gamma_xe*Sigma_f*flux;
		b = gamma_I*Sigma_f*flux;
		c = lambda_xe + sigma_a*flux;
    	N_xe = a/c + b/c - b/(c - lambda_I)*exp(-lambda_I*t) - 
			(a/(c - lambda_I) + (a + b)/c)*exp(-c*t);
    	N_I = b/lambda_I*(1. - exp(-lambda_I*t));
		
    	std::cout << N_xe << " " << sol(0) << std::endl;
    	std::cout << N_I << " " << sol(1) << std::endl;
    	std::cout << " " << std::endl;

		//assert(isApprox(sol(0), N_xe));
		//assert(isApprox(sol(1), N_I));
		//assert(isApprox(sol(2), 1.0));
	}

}

int main(){

	testSolverTime();
	tankProblem();
	//xenonIodineProblem();
}



