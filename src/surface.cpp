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
// @param molarMass  Molar mass of species [g/mol]
// @param initCon    Initial concentration [kg/m^3]
// @param diffCoeff  Diffusion coefficient [m^2/s]
// @param name      Species name
// @param Transport  bool to set if the speices is to be transport with the
//              velocity field
//*****************************************************************************
void surface::addSpecies(double molarMass, double initCon, double diffCoeff,
  std::string name, bool transport){

  species spec(molarMass, initCon, diffCoeff, name, transport);
  speciesVector.push_back(spec);
}

//*****************************************************************************
// Gets a pointer to the species object on the surface
//
// @param specID  ID of the species
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
