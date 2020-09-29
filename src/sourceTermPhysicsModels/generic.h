//*****************************************************************************
// Author: Zack Taylor 
//
//	Defines the generic source term model
//*****************************************************************************
#include <vector>
#include "physicsModelABC.h"

class generic : public physicsModel {

	// Class attributes
	protected:
	std::vector<double> coeffs;

	public:
	// Constructor
	generic();
	// Sets the coefficients for generic source term.
	void setModel(std::vector<double>);
	// Gets the mass transfer coefficient
	double getMassTransferCoeff(int, scalarData*);

};
