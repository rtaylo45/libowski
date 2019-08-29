//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include "meshCellData.h"
#include <assert.h>

//*****************************************************************************
// Adds a species to the cell
//
// @param molarMass	Molar mass of species [lbm/mol]
// @param initCon		Initial concentration [lbm/ft^3]
//*****************************************************************************
void meshCell::addSpecies(double molarMass, double initCon){

	species spec(molarMass, initCon);
	speciesVector.push_back(spec);

}

//*****************************************************************************
// Gets a pointer to the species object in the cell
//
// @param specID	ID of the species [lbm/ft^3]
//*****************************************************************************
species* meshCell::getSpecies(int specID){
	// Checks to make sure the specID is not out of range
	assert(specID <= speciesVector.size() and specID>= 0);
	return &speciesVector[specID];

}

//*****************************************************************************
// Cleans the species in the cell
//*****************************************************************************
void meshCell::cleanSpecies(){
	speciesVector.clear();
	//for(int i = 0; i < speciesVector.size(); i++){
	//		species spec = speciesVector[i];
	//		spec.clean();
	//}
}
