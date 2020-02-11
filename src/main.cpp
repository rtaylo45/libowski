#include <Eigen/Sparse>
#include <iostream>
#include "matrixExponential.h"
#include "mpiProcess.h"
#include <chrono>
#include <math.h>

using namespace std::chrono;
using namespace Eigen;

void tankProblem(){

	typedef Eigen::Triplet<double> T;
   double x1_0 = 1000.0, x2_0 = 0.0, x3_0 = 0.0;
   double t = 0.0; 
	double x1, x2, x3;
	int steps = 1;
	int myid;
	double totalTime = 50.0;
	double dt = totalTime/steps;
   SparseMatrix<double> A(3,3);
   SparseVector<double> b(3);
	MatrixXd sol;
	std::vector<T> tripletList;
	tripletList.reserve(5);
   matrixExponential *expSolver;

	// Set the A matrix and inition condition vector b
	b.insert(0,0) = x1_0;
	tripletList.push_back(T(0,0,-0.5)); tripletList.push_back(T(1,0,0.5));
	tripletList.push_back(T(1,1,-0.25)); tripletList.push_back(T(2,1,0.25));
	tripletList.push_back(T(2,2,-1./6.)); 
	A.setFromTriplets(tripletList.begin(), tripletList.end());

	// Sets the solver
	expSolver = matrixExponentialFactory::getExpSolver("CRAM");
	myid = mpi.rank;

	for (int i = 0; i < steps; i++){

		t = t + dt;

		sol = expSolver->apply(A, b, t);
    	x1 = x1_0*exp(-t/2.);
    	x2 = -2.*x1_0*exp(-t/2.) + (x2_0 + 2.*x1_0)*exp(-t/4.);
    	x3 = (3./2.)*x1_0*exp(-t/2.) - 3.*(x3_0 + 2.*x1_0)*exp(-t/4) +
    	    (x3_0 - (3./2.)*x1_0 + 3.*(x2_0 + 2.*x1_0))*exp(-t/6.);
	
		if (myid == 0){
			std::cout << x1 << " " << sol(0) << std::endl;
    		std::cout << x2 << " " << sol(1) << std::endl;
    		std::cout << x3 << " " << sol(2) << std::endl;
    		std::cout << " " << std::endl;
		}

	}
	mpi.finalize();
	
}


int main(){

	tankProblem();

}
