#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Sparse>
#include <iostream>
#include <chrono> 

using namespace std::chrono; 
using namespace Eigen;

int main()
{
	for (int i = 0; i < 20; i++){
		int n = 10;
		//MatrixXf A = MatrixXf::Random(n,n);
		SparseMatrix<double> sm1(n,n);
		auto start = high_resolution_clock::now();
		//MatrixXf matEx = A.exp();
		auto end = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(end - start);
		//std::cout << "n: " << n << " Time: " << duration.count() << std::endl;
		std::cout << sm1 << std::endl;
	}
}
