//*****************************************************************************
// Author: Zack Taylor 
//
//	Defines the generic source term model
//*****************************************************************************
#include "generic.h"

//**************************************************************************
// Constructor
//
//**************************************************************************
generic::generic(){
}

//*****************************************************************************
// Sets the coefficients for the model
//
// @param coeffs	Vector of coefficients for the source term [1/s]
//*****************************************************************************
void generic::setModel(std::vector<double> coefficientVector){
	coeffs = coefficientVector;
}

//*****************************************************************************
// Gets the coefficient for a species relation
//
// @param otherSpecID		specID of the mass transfer coefficient to find.
// @param scalarVariables	Pointer to scalar data object for the mesh cell
//*****************************************************************************
double generic::getMassTransferCoeff(int otherSpecID, scalarData* scalarVariables){
	return coeffs[otherSpecID];
}
