//*****************************************************************************
// Author: Zack Taylor
//
// Base utility functions that are used in libowski. These are general math
// and helper functions.
//*****************************************************************************
#ifndef UTILBASE_H
#define UTILBASE_H
#include <cmath>
#include "matrixTypes.h"
#include "vectorTypes.h"

using namespace Eigen;

//*****************************************************************************
// Is approx function to test if two numbers are close to one another
//*****************************************************************************
template <typename derived>
bool isApprox(derived, derived, derived = 1e-5, derived = 1e-8);

//*****************************************************************************
// Binomial Coefficient
//*****************************************************************************
int binomialCoeff(int, int);

//*****************************************************************************
// Factorial
//*****************************************************************************
int factorial(int);

//*****************************************************************************
// Creates a linear line space 
//*****************************************************************************
template <typename T>
std::vector<T> lineSpace(T, T, std::size_t);

#endif
