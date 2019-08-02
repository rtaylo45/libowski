#include "coreSolver.h"
#include <iostream>

using namespace Eigen;
//Eigen::Matrix<std::complex<double>,8,1> theta;

// Solver
MatrixXd SolverType::solve(MatrixXd A, MatrixXd w0, double t){

	int s = 8;
	auto At = A*t; // using auto cause i don't know what type this is
	MatrixXd w = 0*w0;
	MatrixXcd ident = MatrixXcd::Identity(A.rows(), A.cols());

	for (int k = 0; k < s; k++){
		auto tempA = At - theta(k)*ident;
		//w = w + 
		//std::cout << theta(k)*ident << std::endl;
		//std::cout << "" << std::endl;
	}
	
	return A;

}
