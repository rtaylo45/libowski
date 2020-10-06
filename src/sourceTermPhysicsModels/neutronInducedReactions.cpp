//*****************************************************************************
// Author: Zack Taylor 
//
//	Defines the source term for neutron induced reactions
//*****************************************************************************
#include "neutronInducedReactions.h"

//**************************************************************************
// Constructor
//
//**************************************************************************
neutronInducedReactions::neutronInducedReactions(){
}

//*****************************************************************************
// Sets the coefficients for the model
//
// @param coefficientVector	Vector of cross sections [cm^2]
//*****************************************************************************
void neutronInducedReactions::setModel(std::vector<double> coefficientVector){
	crossSections = coefficientVector;
}

//*****************************************************************************
// Gets the coefficient for a species relation
//
// @param otherSpecID		specID of the mass transfer coefficient to find.
// @param scalarVariables	Pointer to scalar data object for the mesh cell
//*****************************************************************************
double neutronInducedReactions::getTransitionCoeff(int otherSpecID, 
	scalarData* scalarVariables){
	return crossSections[otherSpecID]*scalarVariables->getNeutronFlux();
}
