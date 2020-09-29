//*****************************************************************************
// Author: Zack Taylor 
//
//*****************************************************************************
#ifndef PHYSICSMODELABC_H 
#define PHYSICSMODELABC_H
#include "scalarData.h"

//*****************************************************************************
// Abstract base for mass transfer models
//*****************************************************************************
#include <vector>
#include "scalarData.h"

class physicsModel {

	// Class methods
	public:
	// Constructor
	physicsModel();	
	// Sets parameters for the generic and transmutation physics models
	virtual void setModel(std::vector<double>) = 0;
	// gets the coefficient
	virtual double getMassTransferCoeff(int, scalarData*) = 0;
	// Sets coefficients for the generic and transmutation physics models
	//virtual void setModel(std::vector<int>, std::vector<double>) = 0;


};

#endif
