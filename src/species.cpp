#include "species.h"
//**************************************************************************
// Constructor
//
// @param molarMass	Molar mass of species [g/mol]
// @param initCon		Initial concentration [kg/m^3]
// @param diffCoeff	Diffusion coefficient [m^2/s]
// @param name			String of the species name
// @param transport	Bool to set if the species is to be transported by the 
//							velocity field
//**************************************************************************
species::species(double molarMass, double initCon, double diffCoeff, 
	std::string name_, bool transport_){
	MM = molarMass;
	c = initCon;
	D = diffCoeff;
	name = name_;
	transport = transport_;
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
	transCoeffs.clear();
}
