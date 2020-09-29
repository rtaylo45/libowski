//*****************************************************************************
// Author: Zack Taylor 
//
//	Defines the source term for neutron induced reactions
//*****************************************************************************
#ifndef PHYSICSMODELNIR_H
#define PHYSICSMODELNIR_H
#include <vector>
#include "physicsModelABC.h"

class neutronInducedReactions : public physicsModel {

	// Class attributes
	protected:
	// Vector of reaction cross section [cm^2]
	std::vector<double> crossSections;

	public:
	// Constructor
	neutronInducedReactions();
	// Sets the coefficients for reaction cross sections
	void setModel(std::vector<double>);
	// Gets the mass transfer coefficient
	double getMassTransferCoeff(int, scalarData*);
};

#endif
