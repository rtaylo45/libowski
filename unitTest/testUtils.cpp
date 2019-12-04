#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include <iostream>
#include <vector>
#include "linearAlgebra.h"
#include "mpiProcess.h"
#include "matrixTypes.h"
#include "vectorTypes.h"

//*****************************************************************************
// Builds a nonsymmetric square matrix sparse
//
// @param n		Matrix size
//*****************************************************************************
SparseMatrixD buildAMatrixSparse(int n){
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(3*n);
   SparseMatrixD A(n,n);

	for (int j = 0; j < n-1; j++){
		tripletList.push_back(T(j,j,-3.0));
	}
	for (int j = 0; j < n-2; j++){
		tripletList.push_back(T(j,j+1,1.0));
	}
	for (int j = n; j < n-1; j++){
		tripletList.push_back(T(j,j-1,1.0));
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	return A;
}
//*****************************************************************************
// Builds a tridiagonal square matrix sparse
//
// @param n		Matrix size
//*****************************************************************************
SparseMatrixD buildJMatrixSparse(int n){
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(3*n);
   SparseMatrixD A(n,n);
	tripletList.push_back(T(0,0,1.0));
	tripletList.push_back(T(0,1,-1.0));

	for (int j = 1; j < n-1; j++){
		tripletList.push_back(T(j,j,2.0/(pow(n,2.0))));
	}
	for (int j = 1; j < n-1; j++){
		tripletList.push_back(T(j,j+1,-1.0));
		tripletList.push_back(T(j,j-1,-1.0));
	}
	tripletList.push_back(T(n-1,n-1,1.0));
	tripletList.push_back(T(n-1,n-2,-1.0));

	A.setFromTriplets(tripletList.begin(), tripletList.end());
	return A;
}

//*****************************************************************************
// Builds a nonsymmetric square matrix Dense
//
// @param n		Matrix size
//*****************************************************************************
MatrixD buildAMatrixDense(int n){
	MatrixD A(n,n);

	for (int j = 0; j < n-1; j++){
		A(j,j) = 3.0;
	}
	for (int j = 0; j < n-2; j++){
		A(j,j+1) = 1.0; 
	}
	for (int j = n; j < n-1; j++){
		A(j,j-1) = 1.0;
	}
	return A;
}
//*****************************************************************************
// Builds a tridiagonal square matrix dense
//
// @param n		Matrix size
//*****************************************************************************
MatrixD buildJMatrixDense(int n){
   MatrixD A(n,n);
	A(0,0) = 1.0;
	A(0,1)= -1.0;

	for (int j = 1; j < n-1; j++){
		A(j,j) = 2.0/(pow(n,2.0));
	}
	for (int j = 1; j < n-1; j++){
		A(j,j+1) = -1.0;
		A(j,j-1) = -1.0;
	}
	A(n-1,n-1) = 1.0;
	A(n-1,n-2) = -1.0;

	return A;
}

//*****************************************************************************
// Sets the pseudo inverse
//
//*****************************************************************************
void testPseudoInverse(int myid){

	typedef Eigen::Triplet<double> T;
	double lambda_I = 2.11E-5;
	double lambda_xe = 2.9306E-5;
	double sigma_a = 2.002E-22;
	double Sigma_f = 9.7532E-1;
	double flux = 2.5E16;
	double gamma_xe = 0.002468;
	double gamma_I = 0.063033;
   SparseMatrixD AxeISparse(3,3);
   SparseMatrixD AtankSparse(3,3);
   SparseMatrixD AaMatrixSparse(5,5);
   SparseMatrixD AjMatrixSparse(5,5);
	MatrixD AtankDense(3,3);
	MatrixD AxeIDense(3,3);
	MatrixD AaMatrixDense(5,5);
	MatrixD AjMatrixDense(5,5);
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

	// Dense matrix for xenon problem
	AxeIDense(0,0) = -lambda_xe - sigma_a*flux; AxeIDense(0,1) = lambda_I; 
	AxeIDense(0,2) = gamma_xe*Sigma_f*flux; AxeIDense(1,1) = -lambda_I; 
	AxeIDense(1,2) = gamma_I*Sigma_f*flux;

	// Builds the dense and sparse matrix 
	AaMatrixDense = buildAMatrixDense(AaMatrixDense.rows());
	AjMatrixDense = buildJMatrixDense(AjMatrixDense.rows());
	AaMatrixSparse = buildAMatrixSparse(AaMatrixSparse.rows());
	AjMatrixSparse = buildJMatrixSparse(AjMatrixSparse.rows());

	std::cout << "Tank Matrix"<<std::endl;
	std::cout << "Dense Matrix A-1*A" << std::endl;
	std::cout << AtankDense.inverse()*AtankDense << std::endl;
	std::cout << "Sparse Matrix A-1*A" <<std::endl;
	std::cout << MoorePenroseInv(AtankSparse)*AtankSparse << std::endl;

	std::cout << "Xenon I Matrix"<<std::endl;
	std::cout << "Dense Matrix A-1*A" <<std::endl;
	std::cout << AxeIDense.completeOrthogonalDecomposition().pseudoInverse()*AxeIDense << std::endl;
	std::cout << "Sparse Matrix A-1*A" <<std::endl;
	std::cout << MoorePenroseInv(AxeISparse)*AxeISparse << std::endl;

	std::cout << "A Matrix"<<std::endl;
	std::cout << "Dense Matrix A-1*A" <<std::endl;
	std::cout << AaMatrixDense.completeOrthogonalDecomposition().pseudoInverse()*AaMatrixDense << std::endl;
	std::cout << "Sparse Matrix A-1*A" <<std::endl;
	std::cout << MoorePenroseInv(AaMatrixSparse)*AaMatrixSparse << std::endl;

	std::cout << "J Matrix"<<std::endl;
	std::cout << "Dense Matrix A-1*A" <<std::endl;
	std::cout << AjMatrixDense.completeOrthogonalDecomposition().pseudoInverse()*AjMatrixDense << std::endl;
	std::cout << "Sparse Matrix A-1*A" <<std::endl;
	std::cout << MoorePenroseInv(AjMatrixSparse)*AjMatrixSparse << std::endl;
}

//*****************************************************************************
// Tests matrix squaring 
//
//*****************************************************************************
void testMatrixSquare(int myid){

	typedef Eigen::Triplet<double> T;
   SparseMatrixD A(2,2);
	MatrixXd Solution(2,2);
	std::vector<T> tripletList;
	tripletList.reserve(4);
	int nMax = 5;

	A.insert(0,0) = 2.;
	A.insert(0,1) = 1.;
	A.insert(1,0) = 1.;
	A.insert(1,1) = 2.;

	for (int n = 0; n < nMax; n++){
		int nPow = pow(2,n);
		// Got from online solution
		// https://yutsumura.com/how-to-find-a-formula-of-the-power-of-a-matrix/
		Solution(0,0) = 0.5*(pow(-1.,nPow) + pow(3,nPow));
		Solution(0,1) = 0.5*(pow(-1.,nPow+1) + pow(3,nPow));
		Solution(1,0) = 0.5*(pow(-1.,nPow+1) + pow(3,nPow));
		Solution(1,1) = 0.5*(pow(-1.,nPow) + pow(3,nPow));

		std::cout << Solution << std::endl;
		std::cout << MatrixSquare(A, nPow) << std::endl;
	}

}
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	//testPseudoInverse(myid);
	testMatrixSquare(myid);

	mpi.finalize();
}
