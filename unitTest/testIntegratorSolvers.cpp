#include <assert.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>

#include "ODEintegrator.h"
#include "mpiProcess.h"
#include "matrixTypes.h"
#include "vectorTypes.h"
#include "utilBase.h"

void tankProblem(int myid, ODEintegrator *intSolver){
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
	int steps = 2000;
	double totalTime = 20.0;
	double maxRelativeError = 0.0;
	double dt = totalTime/steps;
   SparseMatrixD A(3,3);
   VectorD b(3);
	MatrixD sol;
	std::vector<T> tripletList;
	tripletList.reserve(5);

	// Set the A matrix and inition condition vector b
	b(0) = x1_0, b(1) = 0.0, b(2) = 0.0;
	tripletList.push_back(T(0,0,-0.5)); tripletList.push_back(T(1,0,0.5));
	tripletList.push_back(T(1,1,-0.25)); tripletList.push_back(T(2,1,0.25));
	tripletList.push_back(T(2,2,-1./6.)); 
	A.setFromTriplets(tripletList.begin(), tripletList.end());

	for (int i = 0; i < steps; i++){

		t = t + dt;
		sol = intSolver->integrate(A, b, dt);
    	x1 = x1_0*exp(-t/2.);
    	x2 = -2.*x1_0*exp(-t/2.) + (x2_0 + 2.*x1_0)*exp(-t/4.);
    	x3 = (3./2.)*x1_0*exp(-t/2.) - 3.*(x3_0 + 2.*x1_0)*exp(-t/4) +
    	    (x3_0 - (3./2.)*x1_0 + 3.*(x2_0 + 2.*x1_0))*exp(-t/6.);

		b = sol;	
		if (myid==0){
			//std::cout << x1 << " " << sol(0) << std::endl;
    		//std::cout << x2 << " " << sol(1) << std::endl;
    		//std::cout << x3 << " " << sol(2) << std::endl;
    		//std::cout << " " << std::endl;
			//std::cout << abs(x1-sol(0))/x1 << std::endl;
			//std::cout << abs(x2-sol(1))/x2 << std::endl;
			//std::cout << abs(x3-sol(2))/x3 << std::endl;
			maxRelativeError = std::max(abs(x1-sol(0))/x1, maxRelativeError);
			maxRelativeError = std::max(abs(x2-sol(1))/x2, maxRelativeError);
			maxRelativeError = std::max(abs(x3-sol(2))/x3, maxRelativeError);
    		//std::cout << " " << std::endl;

			//assert(isApprox(x1, sol(0), 1.e-10, 1.e-11));
			//assert(isApprox(x2, sol(1), 1.e-10, 1.e-11));
			//assert(isApprox(x3, sol(2), 1.e-10, 1.e-11));
		}
	}
	std::cout << "Tank problem " << maxRelativeError << std::endl;
}

void xenonIodineProblem(int myid, ODEintegrator *intSolver){
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
	double a, b, d, k;
	int steps = 100;
	double totalTime = 10000.0;
	double dt = totalTime/steps;
	double lambda_I = 2.11E-5;
	double lambda_xe = 2.9306E-5;
	double sigma_a = 2.002E-22;
	double Sigma_f = 9.7532E-1;
	double flux = 2.5E16;
	double gamma_xe = 0.002468;
	double gamma_I = 0.063033;
	double xenonMM = 135.0, iodineMM = 135.0;
	double AvogNum = 6.02214076E23;
	double maxRelativeError = 0.0;
   SparseMatrixD A(3,3);
   VectorD N0(3);
	MatrixD sol;
	std::vector<T> tripletList;
	tripletList.reserve(5);

	N0(0) = 0.0, N0(1) = 0.0, N0(2) = N_d_0;
	tripletList.push_back(T(0,0,-lambda_xe - sigma_a*flux));
	tripletList.push_back(T(0,1,lambda_I*xenonMM/iodineMM)); 
	tripletList.push_back(T(0,2,gamma_xe*Sigma_f*flux*xenonMM/AvogNum)); 
	tripletList.push_back(T(1,1,-lambda_I));
	tripletList.push_back(T(1,2,gamma_I*Sigma_f*flux*iodineMM/AvogNum)); 
	A.setFromTriplets(tripletList.begin(), tripletList.end());

	for (int i = 0; i < steps; i++){

		t = t + dt;

		sol = intSolver->integrate(A, N0, dt);

		a = lambda_xe + sigma_a*flux;
		b = gamma_I*Sigma_f*flux*iodineMM/AvogNum;
		d = lambda_I*N_I_0;
		k = N_xe_0 - (d-b)/(a - lambda_I) - 
			(b + gamma_xe*Sigma_f*flux*xenonMM/AvogNum)/a;

		// Xenon solution
    	N_xe = -b/(a-lambda_I)*exp(-lambda_I*t) + b/a +
			d*exp(-lambda_I*t)/(a - lambda_I) + k*exp(-a*t) +
			gamma_xe*Sigma_f*flux*xenonMM/AvogNum/a;

		// Iodine solution
    	N_I = b/lambda_I*(1. - exp(-lambda_I*t)) + N_I_0*exp(-lambda_I*t);

		N0 = sol;
		
		if (myid==0){
			//std::cout << N_xe << " " << sol(0) << std::endl;
    		//std::cout << N_I << " " << sol(1) << std::endl;
			//std::cout << abs(N_xe-sol(0))/N_xe << std::endl;
			//std::cout << abs(N_I-sol(1))/N_I << std::endl;
    		//std::cout << " " << std::endl;
			maxRelativeError = std::max(abs(N_xe-sol(0))/N_xe, maxRelativeError);
			maxRelativeError = std::max(abs(N_I-sol(1))/N_I, maxRelativeError);

		}
	}
	std::cout << "Xenon problem " << maxRelativeError << std::endl;
}

int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;
	ODEintegrator *intSolver;
	std::vector<std::string> solvers {"forward euler", "explicit midpoint", "heun second-order",
		"ralston second-order", "kutta third-order", "heun third-order", "ralston third-order",
		"SSPRK3", "classic fourth-order"};
	std::vector<std::string> methods {"explicit"};

	// Loop over methods
	for (std::string &method : methods){
		// Loops over different solvers
		for (std::string &solverType : solvers){
			std::cout << solverType << std::endl;
			intSolver = integratorFactory::getIntegrator(method, solverType);
			tankProblem(myid, intSolver);
			xenonIodineProblem(myid, intSolver);
			std::cout << " " << std::endl;
		}
	}
	mpi.finalize();
}
