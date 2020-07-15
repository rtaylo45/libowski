//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include "surface.h"

//*****************************************************************************
// Constructure
//*****************************************************************************
surface::surface(){};

//*****************************************************************************
// Initilizes the surface
//*****************************************************************************
void surface::set(){
	isInit = true;
}
//*****************************************************************************
// Adds a species to the surface
//
// @param molarMass	Molar mass of species [lbm/mol]
// @param initCon		Initial concentration [lbm/ft^3]
// @param diffCoeff	Diffusion coefficient [ft^2/s]
// @param name			Species name
//*****************************************************************************
void surface::addSpecies(double molarMass, double initCon, double diffCoeff,
	std::string name){

	species spec(molarMass, initCon, diffCoeff, name);
	speciesVector.push_back(spec);
}

//*****************************************************************************
// Gets a pointer to the species object on the surface
//
// @param specID	ID of the species
//*****************************************************************************
species* surface::getSpeciesPtr(int specID){
	// Checks to make sure the specID is not out of range
	assert(specID <= speciesVector.size() and specID>= 0);
	return &speciesVector[specID];
}

//*****************************************************************************
// Clean
//*****************************************************************************
void surface::clean(){
	isInit = false;
	speciesVector.clear();
}
