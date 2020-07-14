#include "species.h"
//**************************************************************************
// Constructor
//
// @param molarMass	Molar mass of species [lbm/mol]
// @param initCon		Initial concentration [lbm/ft^3]
// @param diffCoeff	Diffusion coefficient [ft^2/s]
// @param name			String of the species name
//**************************************************************************
species::species(double molarMass, double initCon, double diffCoeff, std::string name_){
	MM = molarMass;
	c = initCon;
	D = diffCoeff;
	name = name_;
}

//**************************************************************************
// Clean
//**************************************************************************
void species::clean(){
	c = 0.0;
	MM = 0.0;
	s = 0.0;
	D = 0.0;
	name = "None";
	coeffs.clear();
}
