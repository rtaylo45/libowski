//*****************************************************************************
// Author: Zack Taylor
//
// Base utility functions that are used in libowski. These are general math
// and helper functions.
//*****************************************************************************
#include "utilBase.h"
#include <iostream>

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

	double diff = abs(goalVal - testVal);
	if (diff < rtol) { rtolBool = true; }
	if (diff/goalVal < atol) { atolBool = true; }
	if (rtolBool and atolBool) { retBool = true; }
	return retBool;
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
std::vector<int> lineSpace(int start, int end, std::size_t N){
	int h = (end - start) / static_cast<int>(N-1);
   std::vector<int> xs(N);
   std::vector<int>::iterator x;
   int val;
   for (x = xs.begin(), val = start; x != xs.end(); ++x, val += h) {
       *x = val;
   }
   return xs;
}

// Data types that can use the template functions
template bool isApprox(double goalVal, double testVal, double rtol, double atol);
template bool isApprox(float goalVal, float testVal, float rtol, float atol);
