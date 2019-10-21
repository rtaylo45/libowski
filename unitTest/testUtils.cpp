#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include <iostream>
#include <vector>
#include "linearAlgebra.h"
#include "mpiProcess.h"

void testPseudoInverse(int myid){

	typedef Eigen::Triplet<double> T;
	double lambda_I = 2.11E-5;
	double lambda_xe = 2.9306E-5;
	double sigma_a = 2.002E-22;
	double Sigma_f = 9.7532E-1;
	double flux = 2.5E16;
	double gamma_xe = 0.002468;
	double gamma_I = 0.063033;
   SparseMatrix<double> AxeISparse(3,3);
   SparseMatrix<double> AtankSparse(3,3);
   SparseMatrix<double> AtestSparse(2,3);
	MatrixXd AtankDense(3,3);
	MatrixXd AxeIDense(3,3);
	std::vector<T> tripletListTank;
	std::vector<T> tripletListXeI;
	std::vector<T> tripletListTest;
	tripletListTank.reserve(5);
	tripletListXeI.reserve(5);
	tripletListTest.reserve(6);

	// Builds the sparse matrix for tank problem
	tripletListTank.push_back(T(0,0,-0.5)); tripletListTank.push_back(T(1,0,0.5));
	tripletListTank.push_back(T(1,1,-0.25)); tripletListTank.push_back(T(2,1,0.25));
	tripletListTank.push_back(T(2,2,-1./6.)); 
	AtankSparse.setFromTriplets(tripletListTank.begin(), tripletListTank.end());

	// Builds the sparse matrix for xe I problem
	tripletListXeI.push_back(T(0,0,-lambda_xe - sigma_a*flux));
	tripletListXeI.push_back(T(0,1,lambda_I)); 
	tripletListXeI.push_back(T(0,2,gamma_xe*Sigma_f*flux)); 
	tripletListXeI.push_back(T(1,1,-lambda_I));
	tripletListXeI.push_back(T(1,2,gamma_I*Sigma_f*flux)); 
	AxeISparse.setFromTriplets(tripletListXeI.begin(), tripletListXeI.end());

	// Dense matrix for tank problem
	AtankDense(0,0) = -0.5; AtankDense(1,0) = 0.5; AtankDense(1,1) = -0.25; 
	AtankDense(2,2) = -1./6.; AtankDense(2,1) = 0.25;

	AxeIDense(0,0) = -lambda_xe - sigma_a*flux; AxeIDense(0,1) = lambda_I; 
	AxeIDense(0,2) = gamma_xe*Sigma_f*flux; AxeIDense(1,1) = -lambda_I; 
	AxeIDense(1,2) = gamma_I*Sigma_f*flux;

	std::cout << AtankDense.inverse() << std::endl;
	std::cout << MoorePenroseInv(AtankSparse) << std::endl;

	std::cout << AxeIDense.completeOrthogonalDecomposition().pseudoInverse() << std::endl;
	std::cout << MoorePenroseInv(AxeISparse) << std::endl;

}

int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	testPseudoInverse(myid);

	mpi.finalize();
}
