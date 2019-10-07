#include "species.h"
//**************************************************************************
// Constructor
//
// @param molarMass	Molar mass of species [lbm/mol]
// @param initCon		Initial concentration [lbm/ft^3]
// @param diffCoeff	Diffusion coefficient [ft^2/s]
//**************************************************************************
species::species(double molarMass, double initCon, double diffCoeff){
	MM = molarMass;
	c = initCon;
	D = diffCoeff;
}

//**************************************************************************
// Clean
//**************************************************************************
void species::clean(){
	c = 0.0;
	MM = 0.0;
	s = 0.0;
	D = 0.0;
	coeffs.clear();
}
