//*****************************************************************************
// Author: Zack Taylor 
//
//	Defines the gas mass transfer souce term.
//
//	Model:
//		m_dot = dC_gas/dt = (h*A/V)*(C_liquid_bulk - C_liquid_eq)
//				= dC_liquid/dt = (h*A/V)*(C_liquid_eq - C_liquid_bulk)
//
//		m_dot				= Mass transfer rate [kg/m^3/s] 
//		h					= Mass transfer coefficient [m/s]
//		A					= Area of the gas phase surface in the cell to which mass 
//							  transfer is taking place [m^2]
//		V					= Volumn of the system [m^3]
//		C_liquid_bulk	= Concentration of the species in the bulk liquid 
//								phase [kg/m^3]
//		C_liquid_eq		= Concentration of the speices in equilbrium with the 
//							  liquid phase [kg/m^3]
//		
//		Using the idea gas law, the C_liquid_eq is calculated using the 
//		dimensionless Henrys constant. C_liquid_eq = H*R*T*C_gas, where C_gas
//		is the cocentration of the speices in the gas phase with relation to 
//		the gas volume and not the cell volume. These then need to be divided
//		by the volume of the phase.
//
//		If the current species ID is for the gas phase, the resulting 
//		coefficients are:
//
//		otherSpecID = gasID:		-(h*A/V)*1/(gas_void_fraction*molar_mass)*1000*H*R*T
//		otherSpecID = liquidID:	(h*A/V)*1/(liq_void_fraction)
//
//		If the current species ID is for the liquid phase, the resulting 
//		coefficients are:
//
//		otherSpecID = gasID:		(h*A/V)*1/(gas_void_fraction*molar_mass)*1000*H*R*T
//		otherSpecID = liquidID:	-(h*A/V)*1/(liq_void_fraction)
//
//
//*****************************************************************************
#ifndef PHYSICSMODELGASSPARGING_H
#define PHYSICSMODELGASSPARGING_H
#include <vector>
#include <assert.h>
#include "physicsModelABC.h"
#include "constants.h"

class gasSparging : public physicsModel {

	// Class attributes
	protected:
	// Mass transfer coefficient [m/s]. Possitive means that the transfer
	// is from the liquid to the gas. Negative is from the gas to the liquid
	double massTransferCoeff;
	// Henrys law constant [mol/m^3/Pa]
	double H;
	// species molar mass [g/mol]
	double mm;
	// species ID of the species object that owns me
	int myID;
	// Parent species ID, this is the one that is the source term
	int liquidID;	
	// Daughter species ID, this would be the wall species ID
	int gasID;

	public:
	// Class methods
	gasSparging();
	// Sets the coefficients and ids for the model
	void setModel(double, double, double, int, int, int);
	// Gets the transition coefficient
	double getTransitionCoeff(int, scalarData*);

};

#endif
