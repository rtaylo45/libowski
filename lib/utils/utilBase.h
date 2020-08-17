//*****************************************************************************
// Author: Zack Taylor
//
// Base utility functions that are used in libowski. These are general math
// and helper functions.
//*****************************************************************************
#ifndef UTILBASE_H
#define UTILBASE_H
#include <cmath>
#include <iostream>
#include <fstream>
#include "sys.h"
#include "matrixTypes.h"
#include "vectorTypes.h"

using namespace Eigen;

//*****************************************************************************
// Is approx function to test if two numbers are close to one another
//*****************************************************************************
template <typename derived>
bool isApprox(derived, derived, derived = 1e-5, derived = 1e-8);

//*****************************************************************************
// Finds and repaces all coefficients of a given value in a dense matrix
//*****************************************************************************
template <typename derived>
void findReplace(Matrix<derived, Dynamic, Dynamic>&, derived, derived);

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

//*****************************************************************************
// Writes matrix to CSV
//*****************************************************************************
template <typename derived>
void writeCSV(const Matrix<derived, Dynamic, Dynamic>&, const std::string);

//*****************************************************************************
// Writes vector to CSV
//*****************************************************************************
template <typename derived>
void writeCSV(const Matrix<derived, Dynamic, 1>&, const std::string);

//*****************************************************************************
// Reads a CSV file into a matrix
//*****************************************************************************
template <typename derived>
void readCSV(Matrix<derived, Dynamic, Dynamic>&, const std::string);

//*****************************************************************************
// Reads a CSV file into a matrix
//*****************************************************************************
template <typename derived>
void readCSV(Matrix<derived, Dynamic, 1>&, const std::string);

//*****************************************************************************
// Computes the relative RMSE between two matrices 
//*****************************************************************************
template <typename derived>
derived computeRelativeRMSE(const Matrix<derived, Dynamic, Dynamic>&, 
	const Matrix<derived, Dynamic, Dynamic>&);

//*****************************************************************************
// Computes the relative RMSE between two vectors
//*****************************************************************************
template <typename derived>
derived computeRelativeRMSE(const Matrix<derived, Dynamic, 1>&, 
	const Matrix<derived, Dynamic, 1>&);

#endif
