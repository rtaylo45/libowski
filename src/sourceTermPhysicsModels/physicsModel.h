//*****************************************************************************
// Author: Zack Taylor 
//
// Gas sparging class that defines volumetric source terms for mass transfer
//*****************************************************************************
#ifndef MASSTRANSFER_H 
#define MASSTRANSFER_H

//*****************************************************************************
// Abstract base for mass transfer models
//*****************************************************************************
class massTransfer {

	// Class methods
	public:
	// Constructor
	massTransfer();	
	// Computes the mass transfer coefficient [1/s]
	virtual getMassTransferCoefficient();
		


};

#endif
