//*****************************************************************************
// Author: Zack Taylor 
//
//*****************************************************************************
#ifndef PHYSICSMODELABC_H 
#define PHYSICSMODELABC_H

//*****************************************************************************
// Abstract base for mass transfer models
//*****************************************************************************
class physicsModels {

	// Class methods
	public:
	// Constructor
	physicsModels();	

	// gets the coefficient
	virtual double getMassTransferCoeff() = 0;

};

#endif
