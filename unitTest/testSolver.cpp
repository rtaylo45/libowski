#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include "CRAM.h"
#include "mpiProcess.h"
#include <chrono>
#include <assert.h>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <math.h>

using namespace std::chrono;
using namespace Eigen;

//*****************************************************************************
// Functions used in unit testing
//*****************************************************************************
//
//*****************************************************************************
// Test if two number are approx equal
//
// @param goalVal		The real solution value
// @param testVal		Test value
// @param rtol			Relative tolerance
// @param atol			Absolution tolerance
//
// Values for rtol and atol were taken from the default values for numpys 
// isApprox function.
//*****************************************************************************
bool isApprox(double goalVal, double testVal, double rtol = 1e-5, double atol = 1e-8){
	bool retBool = false;

	double diff = abs(goalVal - testVal);
	if (diff < rtol) { retBool = true; }
	if (diff/goalVal < atol) { retBool = true; }
	return retBool;
}
//*****************************************************************************
// Builds a nonsymmetric square matrix 
//
// @param n		Matrix size
//*****************************************************************************
SparseMatrix<double> buildAMatrix(int n){
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(3*n);
   SparseMatrix<double> A(n,n);

	for (int j = 0; j < n; j++){
		tripletList.push_back(T(j,j,-3.0));
	}
	for (int j = 0; j < n-1; j++){
		tripletList.push_back(T(j,j+1,1.0));
	}
	for (int j = n; j < n; j++){
		tripletList.push_back(T(j,j-1,1.0));
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	return A;
}
//*****************************************************************************
// Builds a tridiagonal square matrix 
//
// @param n		Matrix size
//*****************************************************************************
SparseMatrix<double> buildJMatrix(int n){
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(3*n);
   SparseMatrix<double> A(n,n);
	tripletList.push_back(T(0,0,1.0));
	tripletList.push_back(T(0,1,-1.0));

	for (int j = 1; j < n-1; j++){
		tripletList.push_back(T(j,j,2.0/(pow(n,2.0))));
	}
	for (int j = 1; j < n-1; j++){
		tripletList.push_back(T(j,j+1,-1.0));
		tripletList.push_back(T(j,j+1,-1.0));
	}
	tripletList.push_back(T(n-1,n-1,1.0));
	tripletList.push_back(T(n-1,n-2,-1.0));

	A.setFromTriplets(tripletList.begin(), tripletList.end());
	return A;
}
//*****************************************************************************
// Builds a tridiagonal square matrix 
//
// @param n		Matrix size
//*****************************************************************************
SparseMatrix<double> buildSMatrix(int n){
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(3*n);
   SparseMatrix<double> A(n,n);

	for (int j = 0; j < n; j++){
		tripletList.push_back(T(j,j,3.0));
	}
	for (int j = 0; j < n-1; j++){
		tripletList.push_back(T(j,j+1,-1.0));
	}
	for (int j = 1; j < n; j++){
		tripletList.push_back(T(j,j-1,-1.0));
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	return A;
}

//*****************************************************************************
// Builds a sparse vector of ones
//
// @param n		Vector size
//*****************************************************************************
SparseVector<double> buildN0Vector(int n){
	SparseVector<double> n0(n);
	for (int j = 0; j < n; j++){
		n0.insert(j) = 1.0;
	}
	return n0;

}
//*****************************************************************************
// Writes the time dependent solution for the neutron precursor problem
//
// @param sol	Soltuion at time t
// @param t		Time over which the solve takes place
//*****************************************************************************
void writePrecursorSolution(MatrixXd sol, double t){
	int specs = 6, cells = 16;
	int s = 0, c = 0;
	std::ofstream outputFile;

	//std::string strT = std::to_string(t);
	//std::string fileName = strNumProcs + "procs.out";

	outputFile.open("precursors.out", std::ios_base::app);
	outputFile << "Time: "+std::to_string(t)+"\n";

	for (int i = 0; i < sol.rows()-1; i++){
		outputFile << c << " " << s << " " << sol(i) << std::endl;

		if (s == 5) {s = 0; c++;}
		else {s++;}
	}


}
//*****************************************************************************
// Builds the species matrix for precursor problem
//
//*****************************************************************************
SparseMatrix<double> BuildSpeciesMatrix(MatrixXd coeff, MatrixXd varCoeff,
	int numOfSpecs, int numOfLvls, double flux){

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	int nonZeros = 2*numOfSpecs + 3*numOfSpecs*(numOfLvls-1);
	tripletList.reserve(nonZeros);
   SparseMatrix<double> A(numOfSpecs*numOfLvls+1,numOfSpecs*numOfLvls+1);
	double val1, val2;
	int s = 0; // species index
	int c = 0; // cell  index

	for (int i = 0; i < numOfSpecs*numOfLvls;){
		if (c != 0 and s < numOfSpecs){tripletList.push_back(T(i, i-numOfSpecs, 
				flux));} 

		if (s < numOfSpecs){
			val1 = varCoeff(s) - flux;
			val2 = coeff(c,s);

			tripletList.push_back(T(i, i, val1));
			tripletList.push_back(T(i, A.cols()-1, val2));
			s++;
			i++;
		}
		else{
			s = 0;
			c++;
		}
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	return A;
}

//*****************************************************************************
// Unit test
//*****************************************************************************

void testSolverTime(int myid, int numprocs){
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
	int iters = 20;
	double simTime = 0.0;
	std::ofstream outputFile;

	std::string strNumProcs = std::to_string(numprocs);
	std::string fileName = strNumProcs + "procs.out";
	outputFile.open(fileName);

	for (int n = 100; n <= 1e7; n = n*10){
		if (myid==0){std::cout << "Size: "+std::to_string(n) << std::endl;};
		if (myid==0){outputFile << "Size: "+std::to_string(n)+"\n";};

   	SparseMatrix<double> A(n,n);
   	SparseVector<double> n0(n);
		A = buildAMatrix(n);
		n0 = buildN0Vector(n);

		for (int i = 0; i < iters; i++){
			if (myid==0){std::cout << "Iter: "+std::to_string(i) << std::endl;};

			// Get the solver object
			SolverType ExpSolver;

			// Build random dense matrix
			//MatrixXd A = MatrixXd::Random(n,n);
			//MatrixXd n0 = MatrixXd::Random(n,1);

			MatrixXd sol;

			// Start time
			auto start = high_resolution_clock::now();
			sol = ExpSolver.solve(A, n0, 1.0);
			auto end = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(end - start);

			// Convert to seconds
			if (myid==0){	
				//assert(duration.count()/1.e6 < 0.8);
				simTime = simTime + duration.count()/1.e6;
				outputFile << "Time: "+std::to_string(duration.count()/1.e6)+"\n";
			}

		}
		if (myid==0){
			outputFile << "\n";
		}
	}
	outputFile.close();
}

void tankProblem(int myid){
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
	double totalTime = 20.0;
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
	
		if (myid==0){
			//std::cout << x1 << " " << sol(0) << std::endl;
    		//std::cout << x2 << " " << sol(1) << std::endl;
    		//std::cout << x3 << " " << sol(2) << std::endl;
    		//std::cout << " " << std::endl;
			//std::cout << abs(x1-sol(0))/sol(0) << std::endl;
			//std::cout << abs(x2-sol(1))/sol(1) << std::endl;
			//std::cout << abs(x3-sol(2))/sol(2) << std::endl;

			assert(isApprox(x1, sol(0)));
			assert(isApprox(x2, sol(1)));
			assert(isApprox(x3, sol(2)));
		}
	}
}

void xenonIodineProblem(int myid){
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
	int steps = 1000;
	double totalTime = 10000.0;
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

	N0.insert(2,0) = N_d_0;
	tripletList.push_back(T(0,0,-lambda_xe - sigma_a*flux));
	tripletList.push_back(T(0,1,lambda_I)); 
	tripletList.push_back(T(0,2,gamma_xe*Sigma_f*flux)); 
	tripletList.push_back(T(1,1,-lambda_I));
	tripletList.push_back(T(1,2,gamma_I*Sigma_f*flux)); 
	A.setFromTriplets(tripletList.begin(), tripletList.end());

	for (int i = 0; i < steps; i++){

		t = t + dt;

		sol = ExpSolver.solve(A, N0, t);

		a = lambda_xe + sigma_a*flux;
		b = gamma_I*Sigma_f*flux;
		d = lambda_I*N_I_0;
		k = N_xe_0 - (d-b)/(a - lambda_I) - (b + gamma_xe*Sigma_f*flux)/a;

		// Xenon solution
    	N_xe = -b/(a-lambda_I)*exp(-lambda_I*t) + b/a +
			d*exp(-lambda_I*t)/(a - lambda_I) + k*exp(-a*t) +
			gamma_xe*Sigma_f*flux/a;

		// Iodine solution
    	N_I = b/lambda_I*(1. - exp(-lambda_I*t)) + N_I_0*exp(-lambda_I*t);
		
		if (myid==0){
			std::cout << N_xe << " " << sol(0) << std::endl;
    		std::cout << N_I << " " << sol(1) << std::endl;
			std::cout << abs(N_xe-sol(0))/N_xe << std::endl;
			std::cout << abs(N_I-sol(1))/N_xe << std::endl;
    		std::cout << " " << std::endl;

			assert(isApprox(sol(0), N_xe));
			assert(isApprox(sol(1), N_I));
			
		}
	}
}
void neutronPrecursorProblem(int myid){
//*****************************************************************************
//	Problem statement:
//		Neutron precursors problem. This is a reactor with 16 axial levels, the 
//		generation from fission is varied axially so we need to model each species
//		in each cell with their own set of equations. The total size of the 
//		problem is 6 species by 16 axial levels. 
//
//		dC1/dt = gamma_C1*Sigma_f*flux - lambda_C1*C1
//		dC2/dt = gamma_C2*Sigma_f*flux - lambda_C2*C2
//		dC3/dt = gamma_C3*Sigma_f*flux - lambda_C3*C3
//		dC4/dt = gamma_C4*Sigma_f*flux - lambda_C4*C4
//		dC5/dt = gamma_C5*Sigma_f*flux - lambda_C5*C5
//		dC6/dt = gamma_C6*Sigma_f*flux - lambda_C6*C6
//		dCd/dt = 0.0
//
//	Initial conditons:
//		C1_0 = 0.0
//		C2_0 = 0.0
//		C3_0 = 0.0
//		C4_0 = 0.0
//		C5_0 = 0.0
//		C6_0 = 0.0
//		Cd_0 = 1.0
//
//	To add the constant source terms we need to add a dummy species to hold the
//	coefficients. 
//	
//****************************************************************************
	typedef Eigen::Triplet<double> T;
   double t = 0.0, tnew = 0.0; 
	int steps = 140;
	double totalTime = 140.0;
	double dt = totalTime/steps;
	int numOfSpecs = 6;
	int numOfLvls = 16;
	int nonZeros = 2*numOfSpecs + 3*numOfSpecs*(numOfLvls-1);
	double lambdaC1 = 0.0125, lambdaC2 = 0.0318, lambdaC3 = 0.109;
	double lambdaC4 = 0.3170, lambdaC5 = 1.3500, lambdaC6 = 8.640;
	double val1, val2;
   SparseMatrix<double> A(numOfSpecs*numOfLvls+1,numOfSpecs*numOfLvls+1);
   //SparseVector<double> N0(numOfSpecs*numOfLvls+1);
   VectorXd N0(numOfSpecs*numOfLvls+1);
	MatrixXd sol;
	MatrixXd coeff(16,7);
	VectorXd varCoeff(7);
	double flux = 0.2, diff = 0.0;
	std::vector<T> tripletList;
	tripletList.reserve(nonZeros);

	// Sets the coefficients for decay
	varCoeff(0) = -lambdaC1, varCoeff(1) = -lambdaC2, varCoeff(2) = -lambdaC3;
	varCoeff(3) = -lambdaC4, varCoeff(4) = -lambdaC5, varCoeff(5) = -lambdaC6;
	varCoeff(6) = 0.0;
	
	// Coefficients in order of precursor groups for columns
	coeff(0,0)  = 1.9490E-04, coeff(0,1)  = 1.0149E-03, coeff(0,2)  = 9.8409E-04;
	coeff(1,0)  = 3.5140E-04, coeff(1,1)  = 1.8298E-03, coeff(1,2)  = 1.7743E-03;
	coeff(2,0)  = 4.8659E-04, coeff(2,1)  = 2.5338E-03, coeff(2,2)  = 2.4569E-03;
	coeff(3,0)  = 6.0108E-04, coeff(3,1)  = 3.1300E-03, coeff(3,2)  = 3.0350E-03;
	coeff(4,0)  = 6.9421E-04, coeff(4,1)  = 3.6150E-03, coeff(4,2)  = 3.5052E-03;
	coeff(5,0)  = 7.6455E-04, coeff(5,1)  = 3.9812E-03, coeff(5,2)  = 3.8604E-03;
	coeff(6,0)  = 8.1053E-04, coeff(6,1)  = 4.2207E-03, coeff(6,2)  = 4.0926E-03;
	coeff(7,0)  = 8.3098E-04, coeff(7,1)  = 4.3271E-03, coeff(7,2)  = 4.1958E-03;
	coeff(8,0)  = 8.2533E-04, coeff(8,1)  = 4.2977E-03, coeff(8,2)  = 4.1673E-03;
	coeff(9,0)  = 7.9376E-04, coeff(9,1)  = 4.1333E-03, coeff(9,2)  = 4.0079E-03;
	coeff(10,0) = 7.3720E-04, coeff(10,1) = 3.8388E-03, coeff(10,2) = 3.7223E-03;
	coeff(11,0) = 6.5731E-04, coeff(11,1) = 3.4228E-03, coeff(11,2) = 3.3189E-03;
	coeff(12,0) = 5.5629E-04, coeff(12,1) = 2.8967E-03, coeff(12,2) = 2.8088E-03;
	coeff(13,0) = 4.3661E-04, coeff(13,1) = 2.2735E-03, coeff(13,2) = 2.2045E-03;
	coeff(14,0) = 3.0078E-04, coeff(14,1) = 1.5662E-03, coeff(14,2) = 1.5187E-03;
	coeff(15,0) = 1.5123E-04, coeff(15,1) = 7.8752E-04, coeff(15,2) = 7.6362E-04;

	coeff(0,3)  = 2.8054E-03, coeff(0,4)  = 8.1466E-04, coeff(0,5)  = 2.8792E-04;
	coeff(1,3)  = 5.0582E-03, coeff(1,4)  = 1.4688E-03, coeff(1,5)  = 5.1911E-04;
	coeff(2,3)  = 7.0042E-03, coeff(2,4)  = 2.0339E-03, coeff(2,5)  = 7.1883E-04;
	coeff(3,3)  = 8.6521E-03, coeff(3,4)  = 2.5124E-03, coeff(3,5)  = 8.8795E-04;
	coeff(4,3)  = 9.9927E-03, coeff(4,4)  = 2.9017E-03, coeff(4,5)  = 1.0255E-03;
	coeff(5,3)  = 1.1005E-02, coeff(5,4)  = 3.1957E-03, coeff(5,5)  = 1.1294E-03;
	coeff(6,3)  = 1.1667E-02, coeff(6,4)  = 3.3880E-03, coeff(6,5)  = 1.1974E-03;
	coeff(7,3)  = 1.1961E-02, coeff(7,4)  = 3.4734E-03, coeff(7,5)  = 1.2276E-03;
	coeff(8,3)  = 1.1880E-02, coeff(8,4)  = 3.4498E-03, coeff(8,5)  = 1.2192E-03;
	coeff(9,3)  = 1.1426E-02, coeff(9,4)  = 3.3178E-03, coeff(9,5)  = 1.1726E-03;
	coeff(10,3) = 1.0611E-02, coeff(10,4) = 3.0814E-03, coeff(10,5) = 1.0890E-03;
	coeff(11,3) = 9.4616E-03, coeff(11,4) = 2.7475E-03, coeff(11,5) = 9.7103E-04;
	coeff(12,3) = 8.0074E-03, coeff(12,4) = 2.3252E-03, coeff(12,5) = 8.2178E-04;
	coeff(13,3) = 6.2846E-03, coeff(13,4) = 1.8250E-03, coeff(13,5) = 6.4498E-04;
	coeff(14,3) = 4.3295E-03, coeff(14,4) = 1.2572E-03, coeff(14,5) = 4.4433E-04;
	coeff(15,3) = 2.1769E-03, coeff(15,4) = 6.3215E-04, coeff(15,5) = 2.2341E-04;

	coeff(0,6)  = 0.0, coeff(1,6)  = 0.0, coeff(2,6)  = 0.0, coeff(3,6)  = 0.0; 
	coeff(4,6)  = 0.0, coeff(5,6)  = 0.0, coeff(6,6)  = 0.0, coeff(7,6)  = 0.0; 
	coeff(8,6)  = 0.0, coeff(9,6)  = 0.0, coeff(10,6) = 0.0, coeff(11,6) = 0.0; 
	coeff(12,6) = 0.0, coeff(13,6) = 0.0, coeff(14,6) = 0.0, coeff(15,6) = 0.0; 

	A = BuildSpeciesMatrix(coeff, varCoeff, numOfSpecs, numOfLvls, flux);
	N0(A.cols()-1) = 1.0;
	

	// Sets the solver
   SolverType ExpSolver;

	//auto start = high_resolution_clock::now();
	for (int i = 0; i < 80; i++){

		t = t + dt;

		sol = ExpSolver.solve(A, N0, t);
		if (myid==0){ writePrecursorSolution(sol, t);}
	}

	flux = 0.08;
	N0 = sol;
	A = BuildSpeciesMatrix(0.8*coeff, varCoeff, numOfSpecs, numOfLvls, 
		flux);

	for (int i = 80; i < steps; i++){

		t = t + dt;
		tnew = tnew + dt;

		sol = ExpSolver.solve(A, N0, tnew);
		if (myid==0){ writePrecursorSolution(sol, t);}
	}
	//auto end = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(end - start);
	//std::cout << myid << " " << duration.count()/1.e6 << std::endl;

}

int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;


	//testSolverTime(myid, numprocs);
	tankProblem(myid);
	xenonIodineProblem(myid);
	neutronPrecursorProblem(myid);

	mpi.finalize();
}



