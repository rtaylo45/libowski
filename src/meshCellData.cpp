//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include "meshCellData.h"

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
// Cleans the species in the cell
//*****************************************************************************
void meshCell::cleanSpecies(){
	for(int i = 0; i < speciesVector.size(); i++){
			species spec = speciesVector[i];
			spec.clean();
	}
}
