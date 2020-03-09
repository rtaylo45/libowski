//*****************************************************************************
// Author: Zack Taylor
//
// Functions used for assertions in libowski
//*****************************************************************************
#ifndef UTILASSERTIONS_H
#define UTILASSERTIONS_H
#include "matrixTypes.h"
#include "vectorTypes.h"

using namespace Eigen;


//*****************************************************************************
// Pseudo inverse linearly independent columns for sparse matrix
//*****************************************************************************
template <typename derived>
bool isApprox(derived, derived, derived = 1e-5, derived = 1e-8);
#endif
