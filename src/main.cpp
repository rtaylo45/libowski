//#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Core>
#include <iostream>
#include "coreSolver.h"
#include <chrono>

using namespace std::chrono;
using namespace Eigen;

int main()
{
	
	for (int i = 0; i < 20; i++){
		int n = 2000;
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
