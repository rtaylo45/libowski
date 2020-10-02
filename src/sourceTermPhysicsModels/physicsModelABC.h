//*****************************************************************************
// Author: Zack Taylor 
//
//*****************************************************************************
#ifndef PHYSICSMODELABC_H 
#define PHYSICSMODELABC_H
#include <vector>
#include "scalarData.h"

//*****************************************************************************
// Abstract base for mass transfer models
//*****************************************************************************
class physicsModel {

	// Class methods
	public:
	// Constructor
	physicsModel();	
	// Sets parameters for the generic and transmutation physics models
	virtual void setModel(std::vector<double>){};
	// gets the coefficient
	virtual double getTransitionCoeff(int, scalarData*) = 0;
	// Sets coefficients for the generic and transmutation physics models
	virtual void setModel(std::vector<int>, std::vector<double>){};
	// Sets coefficients for the wall deposition model
	virtual void setModel(double, int, int, int, bool){};

};

#endif
