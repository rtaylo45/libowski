//#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Sparse>
#include <iostream>
#include "coreSolver.h"
#include <chrono>
#include <math.h>

using namespace std::chrono;
using namespace Eigen;

void tankProblem(){

	SolverType ExpSolver;
	MatrixXd A = MatrixXd::Zero(3,3);
	MatrixXd b = MatrixXd::Zero(3,1);
	double x1_0 = 100.0; double x2_0 = 0.0; double x3_0 = 0.0;
	double t = 10.0;

	A(0,0) = -0.5; A(1,0) = 0.5; A(1,1) = -0.25;
	A(2,1) = 0.25; A(2,2) = -1./6.; b(0,0) = x1_0;

	MatrixXd sol = ExpSolver.solve(A, b, 10.0);
	std::cout << A << std::endl;
	std::cout << "" << std::endl;
	std::cout << b << std::endl;
	std::cout << "" << std::endl;
	std::cout << sol << std::endl;
	double x1 = x1_0*exp(-t/2.);
	double x2 = -2.*x1_0*exp(-t/2.) + (x2_0 + 2.*x1_0)*exp(-t/4.);
	double x3 = (3./2.)*x1_0*exp(-t/2.) - 3.*(x3_0 + 2.*x1_0)*exp(-t/4) + 
		(x3_0 - (3./2.)*x1_0 + 3.*(x2_0 + 2.*x1_0))*exp(-t/6.);

	std::cout << x1 << " " << sol(0) << std::endl;
	std::cout << x2 << " " << sol(1) << std::endl;
	std::cout << x3 << " " << sol(2) << std::endl;
}

void loopTimeProblem(){

	for (int i = 0; i < 10; i++){
		int n = 500;
		SolverType ExpSolver;
		MatrixXd A = MatrixXd::Random(n,n);
		MatrixXd w0 = MatrixXd::Random(n,1);

		//std::cout << "A matrix" << std::endl;
		//std::cout << A << std::endl;
		//std::cout << "" << std::endl;
		//std::cout << "w0 vector" << std::endl;
		//std::cout << w0 << std::endl;
		//std::cout << "" << std::endl;

		auto start = high_resolution_clock::now();
		MatrixXd sol = ExpSolver.solve(A, w0, 0.1);
		auto end = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(end - start);
		std::cout << "n: " << n << " Time: " << duration.count()/1.e6 << std::endl;
		//std::cout << "sol" << std::endl;
		//std::cout << sol << std::endl;
		//std::cout << "" << std::endl;
	}
}

int main(){

	//loopTimeProblem();	
	tankProblem();

}
