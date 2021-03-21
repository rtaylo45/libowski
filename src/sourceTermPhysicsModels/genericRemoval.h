//*****************************************************************************
// Author: Zack Taylor 
//
//	Defines the generic removal source term model
//*****************************************************************************
#ifndef PHYSICSMODELGENERICREMOVAL_H
#define PHYSICSMODELGENERICREMOVAL_H
#include "physicsModelABC.h"

class genericRemoval : public physicsModel {

	// Class attributes
	protected:
	// Coefficient for the removal rate [1/s]
	double removalCoeff;
	// The id of the species that needs to removed i.e. the species ID that
	// owns the model
	int myID;

	public:
	// Constructor
	genericRemoval();
	// Sets removal coefficient in units of [1/s]
	void setModel(int, double);
	// Gets the transition coefficient
	double getTransitionCoeff(int, scalarData*);

};

#endif
