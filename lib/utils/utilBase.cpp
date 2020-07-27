//*****************************************************************************
// Author: Zack Taylor
//
// Base utility functions that are used in libowski. These are general math
// and helper functions.
//*****************************************************************************
#include "utilBase.h"

using namespace Eigen;

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
template <typename derived>
bool isApprox(derived goalVal, derived testVal, derived rtol, derived atol){
	bool retBool = false;
	bool rtolBool = false;
	bool atolBool = false;

	double diff = std::abs(goalVal - testVal);
	if (diff < rtol) { rtolBool = true; }
	if (diff/goalVal < atol) { atolBool = true; }
	if (rtolBool and atolBool) { retBool = true; }
	return retBool;
}

//*****************************************************************************
// Find and replaces all coefficients of a given value in a dense matrix.
//
// Need to find the more efficient way to do it
//
// Right now im only using this for an integer matrix so i don't use a tol for
// how close the two numbers are.
//
// @param A					Dense matrix with the values 
// @param findValue		Value to find
// @param replaceValue	Value to replace
//*****************************************************************************
template <typename derived>
void findReplace(Matrix<derived, Dynamic, Dynamic>& A, derived findValue, 
	derived replaceValue){

	for (int i = 0; i < A.rows(); i++){
		for (int j = 0; j < A.cols(); j++){
			if (A(i,j) == findValue){
				A(i,j) = replaceValue;
			}				
		}
	}
}

//*****************************************************************************
// Binomial coefficient. N choose k
// Taken from:
// https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
//
// @param N		Number of things
// @param k		Number to be taken
//*****************************************************************************
int binomialCoeff(int n, int k){  
    int res = 1;  
  
   // Since C(n, k) = C(n, n-k)  
   if ( k > n - k ){
      k = n - k;  
	} 
   // Calculate value of  
   // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]  
   for (int i = 0; i < k; ++i){  
      res *= (n - i);  
      res /= (i + 1);  
   }  
  
   return res;  
} 
//*****************************************************************************
// Factorial
//
// @param n		Factorial number
//*****************************************************************************
int factorial(int n)
{
   if(n > 1){
		return n * factorial(n - 1);
	}
   else{
		return 1;
	}		
}
//*****************************************************************************
// Creates a linear line space 
//
// @param a		Starting point
// @param b		Ending point
// @param N		Number of points
//*****************************************************************************
template <typename T>
std::vector<T> lineSpace(T start, T end, std::size_t N){
	T h = (end - start) / static_cast<T>(N-1);
   std::vector<T> xs(N);
   typename std::vector<T>::iterator x;
   T val;
   for (x = xs.begin(), val = start; x != xs.end(); ++x, val += h) {
       *x = val;
   }
   return xs;
}
//*****************************************************************************
// Writes a matrix to a CSV file
//
//	@param A			Matrix
// @param fname	file name
//*****************************************************************************
template <typename derived>
void writeCSV(const Matrix<derived, Dynamic, Dynamic>& A, const std::string fname){
	const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

	std::ofstream outputFile;
	outputFile.open(fname);
	outputFile << A.format(CSVFormat) << std::endl;

}

//*****************************************************************************
// Reads in a matrix to a CSV file
//
// @param A		Reference to the matrix that hold the read data
// @param path	Path to the csv file
//*****************************************************************************
template <typename derived>
void readCSV(Matrix<derived, Dynamic, Dynamic>& A, const std::string path) {
   std::ifstream indata;
   std::string line;
   std::vector<long double> values;
   int rows = 0, cols = 0;

   indata.open(path);
   while (std::getline(indata, line)) {
      std::stringstream lineStream(line);
      std::string cell;
      while (std::getline(lineStream, cell, ',')) {
			values.push_back(std::stold(cell));
      }
      ++rows;
   }
	cols = values.size()/rows;
	A = Matrix<derived, Dynamic, Dynamic>::Zero(rows, cols);
	// Loops over vector to add it to the matrix
	int tempRow = 0, tempCol = 0;
	for (int i = 0; i < values.size(); i++){
		A(tempRow, tempCol) = (derived)values[i];
		tempCol += 1;
		if (tempCol >= cols){tempCol = 0; tempRow += 1;};
	}
}

//*****************************************************************************
// Computes the relative RMSE between two matrices 
// 
// @param refMat
// @param apporxMat
//*****************************************************************************
template <typename derived>
derived computeRelativeRMSE(const Matrix<derived, Dynamic, Dynamic>& refMat,
	const Matrix<derived, Dynamic, Dynamic>& approxMat){
	// make sure the mats are the same size
	assert(refMat.cols() == approxMat.cols());
	assert(refMat.rows() == approxMat.rows());
	int rows = refMat.rows();
	int cols = refMat.cols();
	int N = rows*cols;
	derived ref, approx, error = 0.0, eps = 1.e-16;

	// Loop over mats	
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			ref = refMat(i,j);
			approx = approxMat(i,j);
			if (ref >= eps){
				error += std::pow((ref - approx)/ref, 2.0);
			}
		}
	}
	error = error/(derived)N;
	error = std::pow(error, 0.5);
	return error;
}


// Data types that can use the template functions
template bool isApprox(double goalVal, double testVal, double rtol, double atol);
template bool isApprox(float goalVal, float testVal, float rtol, float atol);
template std::vector<double> lineSpace(double start, double end, std::size_t N);
template std::vector<int> lineSpace(int start, int end, std::size_t N);
template void findReplace(MatrixLI& A, long int findValue, long int replaceValue);
template void findReplace(MatrixI& A, int findValue, int replaceValue);
template void writeCSV(const MatrixD& A, const std::string fname);
template void writeCSV(const MatrixLD& A, const std::string fname);
template void readCSV(MatrixD& A, const std::string path);
template double computeRelativeRMSE(const MatrixD& refMat, const MatrixD& approxMat);
