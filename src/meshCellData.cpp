//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include "meshCellData.h"

//**************************************************************************
// Constructor
//
// @param iIndex			Index of cell in x direction
// @param jIndex			Index of cell in y direction
// @param absoluteIndex	Absolute index of the cell
// @param xCor				Location of cell center in x direction
// @param yCor				Location of cell center in y direction
// @param dx				dx of cell
// @param dy				dy of cell
//**************************************************************************
meshCell::meshCell(int iIndex, int jIndex, int absoluteIndex, double xCor, 
		double yCor, double dx_, double dy_){
	i = iIndex;
	j = jIndex;
	absIndex = absoluteIndex;
	x = xCor;
	y = yCor;
	dx = dx_;
	dy = dy_;	
}

//*****************************************************************************
// Adds a species to the cell
//
// @param molarMass	Molar mass of species [lbm/mol]
// @param initCon		Initial concentration [lbm/ft^3]
// @param diffCoeff	Diffusion coefficient [ft^2/s]
//*****************************************************************************
void meshCell::addSpecies(double molarMass, double initCon, double diffCoeff){

	species spec(molarMass, initCon, diffCoeff);
	speciesVector.push_back(spec);

}

//*****************************************************************************
// Gets a pointer to the species object in the cell
//
// @param specID	ID of the species
//*****************************************************************************
species* meshCell::getSpecies(int specID){
	// Checks to make sure the specID is not out of range
	assert(specID <= speciesVector.size() and specID>= 0);
	return &speciesVector[specID];

}
//*****************************************************************************
// Gets species concentration
//*****************************************************************************
double meshCell::getSpecCon(int specID){
	assert(specID <= speciesVector.size() and specID>= 0);
	species* spec = getSpecies(specID);
	return spec->c;
}

//*****************************************************************************
// Gets a pointer to the species object in the cell
//
// @param con		Concentration [lbm/ft^3]
// @param specID	ID of the species
//*****************************************************************************
void meshCell::setSpeciesConcentration(double con, int specID){
	species spec = speciesVector[specID];
	spec.c = con;

}

//*****************************************************************************
// Cleans the species in the cell
//*****************************************************************************
void meshCell::cleanSpecies(){
	speciesVector.clear();
}
