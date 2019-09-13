#include "species.h"
//**************************************************************************
// Constructor
//
// @param molarMass	Molar mass of species [lbm/mol]
// @param initCon		Initial concentration [lbm/ft^3]
//**************************************************************************
species::species(double molarMass, double initCon){
	MM = molarMass;
	c = initCon;
}

//**************************************************************************
// Clean
//**************************************************************************
void species::clean(){
	c = 0.0;
	MM = 0.0;
	s = 0.0;
	coeffs.clear();
}
