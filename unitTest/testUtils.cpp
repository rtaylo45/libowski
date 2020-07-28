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
#include "utilBase.h"

//*****************************************************************************
// Prints vector
//*****************************************************************************
void print(std::vector <int> const &a) {

   for(int i=0; i < a.size(); i++)
      std::cout << a.at(i) << ' ';
}

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
		tripletList.push_back(T(j,j,2.0 + 1./(pow(n,2.0))));
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
// Tests the pseudo inverse
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
// Test the factorial function
//*****************************************************************************
void testFactorial(int myid){

	// Test values
	const int testValue1 = 3*2*1;
	const int testValue2 = 4*3*2*1;
	const int testValue3 = 5*4*3*2*1;
	const int testValue4 = 6*5*4*3*2*1;

	// Test the function
	assert(factorial(3) == testValue1);
	assert(factorial(4) == testValue2);
	assert(factorial(5) == testValue3);
	assert(factorial(6) == testValue4);
}

//*****************************************************************************
// Test the binomial coefficient function
//*****************************************************************************
void testBinomialCoeff(int myid){

	// Test values
	const int n1test = 12;
	const int k1test = 2;
	const int test1Answer = 66;

	const int n2test = 50;
	const int k2test = 3;
	const int test2Answer = 19600;

	const int n3test = 12;
	const int k3test = 0;
	const int test3Answer = 1;

	const int n4test = 24;
	const int k4test = 8;
	const int test4Answer = 735471;

	// Test the function
	assert(binomialCoeff(n1test, k1test) == test1Answer);
	assert(binomialCoeff(n2test, k2test) == test2Answer);
	assert(binomialCoeff(n3test, k3test) == test3Answer);
	assert(binomialCoeff(n4test, k4test) == test4Answer);
}

//*****************************************************************************
// Test the Arnoldi algorithm
//
// This unit test exploits the properties:
//
//		AV = VH
//		H = V*HV
//
//	for n = m.
//*****************************************************************************
void testArnoldi(int myid){

	MatrixD V;
	SparseMatrixD H;
	SparseMatrixD A = buildJMatrixSparse(5);
	VectorD b = VectorD::Ones(A.cols());
	double eps = 1.e-12;
	int matrixSize = 5;

	// Test matrix
	arnoldi(A, b, 5, V, H);
	// checks to make sure they are equal
	assert((A*V - V*H).norm() < eps);
	assert((MatrixD(H) - V.adjoint()*A*V).norm() < eps);

	// Loops and builds random matrices
	for (int i = 0; i < 20; i++){
		A = MatrixD::Random(matrixSize,matrixSize).sparseView();
		b = VectorD::Ones(A.cols());
		arnoldi(A, b, matrixSize, V, H);

		// checks to make sure they are equal
		assert((A*V - V*H).norm() < eps);
		assert((MatrixD(H) - V.adjoint()*A*V).norm() < eps);
		
	}
	
}

//*****************************************************************************
// Test the linspace function 
//*****************************************************************************
void testlineSpace(int myid){

	std::vector<int> testVect1, testVect2, testVect3;
	std::vector<int> goalVect1, goalVect2, goalVect3;
	goalVect1 = {1, 2, 3, 4, 5};
	goalVect2 = {2, 4, 6, 8, 10};
	goalVect3 = {20, 40, 60, 80, 100};

	testVect1 = lineSpace(1, 5, 5);
	testVect2 = lineSpace(2, 10, 5);
	testVect3 = lineSpace(20, 100, 5);

	for (int i = 0; i < 5; i++){
		assert(testVect1[i] == goalVect1[i]);
		assert(testVect2[i] == goalVect2[i]);
		assert(testVect3[i] == goalVect3[i]);
	}

}

//*****************************************************************************
// Test the linspace function 
//*****************************************************************************
void testCSVReader(int myid){
	MatrixD testMat1 = MatrixD::Random(4,7);
	MatrixD retMat1;
	double error;

	// writes matrix to csv file
	writeCSV(testMat1, "test1.csv");

	// reads in the csv file
	readCSV(retMat1, "test1.csv");

	// checks error
	error = (testMat1 - retMat1).norm();
	assert(error < 1.e-14);

	std::remove("test1.csv");
}

int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	//testPseudoInverse(myid);
	testFactorial(myid);
	testBinomialCoeff(myid);
	testArnoldi(myid);
	testlineSpace(myid);
	testCSVReader(myid);

	mpi.finalize();
}
