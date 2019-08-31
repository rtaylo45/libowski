//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include "meshCellData.h"
#include <assert.h>

//**************************************************************************
// Constructor
//
// @param iIndex			Index of cell in x direction
// @param jIndex			Index of cell in y direction
// @param absoluteIndex	Absolute index of the cell
// @param xCor				Location of cell center in x direction
// @param yCor				Location of cell center in y direction
//**************************************************************************
meshCell::meshCell(int iIndex, int jIndex, int absoluteIndex, double xCor, 
		double yCor){
	i = iIndex;
	j = jIndex;
	absIndex = absoluteIndex;
	x = xCor;
	y = yCor;	
}

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
}
