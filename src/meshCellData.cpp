//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include "meshCellData.h"

void meshCell::addSpecies(double molarMass, double initCon){

	species spec(molarMass, initCon);
	speciesVector.push_back(spec);

}
