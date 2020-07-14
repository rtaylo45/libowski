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
	// 2D problem
	if (dx_ and dy_){
		volume = dx_*dy_;
	}
	// 1D problem in the x direction
	else if (dx_ and not dy_){
		volume = dx_;
	}
	// 1D problem in the y direction
	else if (dy_ and not dx_){
		volume = dy_;
	}
}

//*****************************************************************************
// Adds a species to the cell
//
// @param molarMass	Molar mass of species [lbm/mol]
// @param initCon		Initial concentration [lbm/ft^3]
// @param diffCoeff	Diffusion coefficient [ft^2/s]
// @param name			Species name
//*****************************************************************************
void meshCell::addSpecies(double molarMass, double initCon, double diffCoeff,
	std::string name){

	species spec(molarMass, initCon, diffCoeff, name);
	// loop over connections to see if the species needs to be added
	// to a surface
	for (int conCount = 0; conCount < connections.size(); conCount ++){
		connection* thisCon = getConnection(conCount);
		surface* conSurface = thisCon->getSurface();
		// if the pointer is not null then add species to that surface
		if(conSurface->isInit){
			conSurface->addSpecies(molarMass, initCon, diffCoeff, name);
		}
	}
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
// Sets the species concentration
//
// @param con		Concentration [lbm/ft^3]
// @param specID	ID of the species
//*****************************************************************************
void meshCell::setSpeciesConcentration(double con, int specID){
	species spec = speciesVector[specID];
	spec.c = con;
}

//*****************************************************************************
// Sets the cells pressure
//
// @param pressure	Pressure in lbf/in^2
//*****************************************************************************
void meshCell::setPressure(double pressure){
	P = pressure;
}

//*****************************************************************************
// Sets the cells temperature
//
// @param temp		Temperature in kelvin
//*****************************************************************************
void meshCell::setTemperature(double temp){
	T = temp;
}

//*****************************************************************************
// Sets the cells scalar neutron flux
//
// @param phi_		Neutron flux in 1/ft^2/s
//*****************************************************************************
void meshCell::setNeutronFlux(double phi_){
	phi = phi_;
}

//*****************************************************************************
// Gets a pointer to a cell connection
//
// @param conID		ID of the connection
//								north = 0
//								south = 1
//								east = 2
//								west = 3
//*****************************************************************************
connection* meshCell::getConnection(int conID){
	// Checks to make sure the conID is not out of range
	assert(conID <= connections.size() and conID>= 0);
	assert(connections[conID].loc == conID);
	return &connections[conID];
}

//*****************************************************************************
// Adds a surface to the cell
//
// @param locID			location ID of the surface to add
//								north = 0
//								south = 1
//								east = 2
//								west = 3
//*****************************************************************************
void meshCell::addSurface(int locID){
	connection* myCon = getConnection(locID);
	myCon->addSurface();
}

//*****************************************************************************
// Cleans the species in the cell
//*****************************************************************************
void meshCell::cleanSpecies(){
	speciesVector.clear();
	// loop over cell connections
	for (int conCount = 0; conCount < connections.size(); conCount ++){
		connection* thisCon = getConnection(conCount);
	}
}
