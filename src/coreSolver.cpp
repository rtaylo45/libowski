#include "coreSolver.h"
#include <iostream>

using namespace Eigen;

// Solver
MatrixXd SolverType::solve(MatrixXd A, MatrixXd w0, double t){

	int s = 8;
	MatrixXd At = A*t;
	MatrixXcd w = 0.*w0;
	MatrixXcd ident = MatrixXcd::Identity(A.rows(), A.cols());

	for (int k = 0; k < s; k++){
		auto tempA = At - theta(k)*ident;
		auto tempB = alpha(k)*w0;

		w = w + tempA.colPivHouseholderQr().solve(tempB);
	}
	w = 2.*w.real();
	w = w + alpha_0*w0;

	return w.real();
}
