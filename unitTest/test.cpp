#include <Eigen/Core>
#include "coreSolver.h"
#include <chrono>
#include <assert.h>
#include <iostream>

using namespace std::chrono;
using namespace Eigen;

bool isApprox(double val1, double val2, double tol = 1e-10){

	double diff = abs(val1 - val2);
	if (diff < tol) { return true; }
	else { return false; }
}

void testSolverTime(){

	for (int i = 1; i < 20; i++){
	int n = 100;

	// Get the solver object
	SolverType ExpSolver;

	// Build random dense matrix
	MatrixXd A = MatrixXd::Random(n,n);
	MatrixXd n0 = MatrixXd::Random(n,1);

	// Start time
	auto start = high_resolution_clock::now();
	MatrixXd sol = ExpSolver.solve(A, n0, 1.0);
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(end - start);

	// Convert to seconds
	assert(duration.count()/1.e6 < 0.5);
	//std::cout << "Time: " << duration.count()/1.e6 << std::endl;

	}
}

void tankProblem(){

    MatrixXd A = MatrixXd::Zero(3,3);
    MatrixXd b = MatrixXd::Zero(3,1);
	MatrixXd sol;
    double x1_0 = 100.0; double x2_0 = 0.0; double x3_0 = 0.0;
    double t = 0.0; double x1; double x2; double x3;
	int steps = 100;
	double totalTime = 100.0;
	double dt = totalTime/steps;

	// Set the A matrix and inition condition vector b
    A(0,0) = -0.5; A(1,0) = 0.5; A(1,1) = -0.25;
    A(2,1) = 0.25; A(2,2) = -1./6.; b(0,0) = x1_0;

	// Sets the solver
    SolverType ExpSolver;

	for (int i = 0; i < steps; i++){

		t = t + dt;

		sol = ExpSolver.solve(A, b, t);
    	x1 = x1_0*exp(-t/2.);
    	x2 = -2.*x1_0*exp(-t/2.) + (x2_0 + 2.*x1_0)*exp(-t/4.);
    	x3 = (3./2.)*x1_0*exp(-t/2.) - 3.*(x3_0 + 2.*x1_0)*exp(-t/4) +
    	    (x3_0 - (3./2.)*x1_0 + 3.*(x2_0 + 2.*x1_0))*exp(-t/6.);
		
    	//std::cout << x1 << " " << sol(0) << std::endl;
    	//std::cout << x2 << " " << sol(1) << std::endl;
    	//std::cout << x3 << " " << sol(2) << std::endl;

		assert(isApprox(sol(0), x1));
		assert(isApprox(sol(1), x2));
		assert(isApprox(sol(2), x3));
	}
}

int main(){

	testSolverTime();
	tankProblem();
}
