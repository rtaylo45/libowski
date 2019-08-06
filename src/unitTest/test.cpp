#include <Eigen/Core>
#include "coreSolver.h"
#include <chrono>
#include <assert.h>
#include <iostream>

using namespace std::chrono;
using namespace Eigen;

void testSolverTime(){

	for (int i = 1; i < 20; i++){
	int n = 1000;

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
	assert(duration.count()/1.e6 < 1.0);
	std::cout << "Time: " << duration.count()/1.e6 << std::endl;

	}
}

int main(){

	testSolverTime();
}
