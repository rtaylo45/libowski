//*****************************************************************************
// Author: Zack Taylor 
//
//	Defines the wall deposition souce term.
//
//	Model:
//		m_dot = dC_surface/dt = (h*A/V)*(C_liquid - C_surface)
//				= dC_liquid/dt = (h*A/V)*(C_surface - C_liquid)
//
//		m_dot	 	 = Mass transfer rate [kg/m^3/s] 
//		h			 = Mass transfer coefficient [m/s]
//		A			 = Area of the surface in the cell to which mass transfer is
//					   taking place [m^2]
//		V			 = Volumn of the system [m^3]
//		C_liquid	 = Concentration of the species in the liquid phase [kg/m^3]
//		C_surface = Concentration of the species on the wall surface [kg/m^3]
//
//	Infinite sink assumption:
//		This assumes that C_surface becomes zero.
//
//		m_dot = (h*A/V)*C_liquid
//*****************************************************************************
#ifndef PHYSICSMODELWALLDEPOSITION_H
#define PHYSICSMODELWALLDEPOSITION_H
#include <vector>
#include "physicsModelABC.h"

class wallDeposition : public physicsModel {

	// Class attributes
	protected:
	// Mass transfer coefficient [m/s]. Possitive 
	// means that the transfer is from the liquid to the wall. Negative
	// is from the wall to the liquid
	double massTransferCoeff;
	// species ID of the species object that owns me
	int myID;
	// Parent species ID, this is the one that is the source term
	int liquidID;	
	// Daughter species ID, this would be the wall species ID
	int surfaceID;
	// Logical to see if an infinite sink assumption is used
	bool infiniteSink = true;

	// Class methods
	public:
	// Constructor
	wallDeposition();
	// Sets the coefficients and ids for the model
	void setModel(double, int, int, int, bool);
	// Gets the transition coefficient
	double getTransitionCoeff(int, scalarData*);

};

#endif
