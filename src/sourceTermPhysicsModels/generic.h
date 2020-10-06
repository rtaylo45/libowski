//*****************************************************************************
// Author: Zack Taylor 
//
//	Defines the generic source term model
//*****************************************************************************
#ifndef PHYSICSMODELGENERIC_H
#define PHYSICSMODELGENERIC_H
#include <vector>
#include "physicsModelABC.h"

class generic : public physicsModel {

	// Class attributes
	protected:
	// Coefficient vector of reaction rates [1/s]
	std::vector<double> coeffs;

	public:
	// Constructor
	generic();
	// Sets the coefficients for generic source term.
	void setModel(std::vector<double>);
	// Gets the transition coefficient
	double getTransitionCoeff(int, scalarData*);

};

#endif
