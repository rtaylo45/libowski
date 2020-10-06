//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include "meshCell.h"

//**************************************************************************
// Constructor
//
// @param iIndex			Index of cell in x direction
// @param jIndex			Index of cell in y direction
// @param absoluteIndex	Absolute index of the cell
// @param xCor				Location of cell center in x direction [m]
// @param yCor				Location of cell center in y direction [m]
// @param dx				dx of cell [m]
// @param dy				dy of cell [m]
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
// @param molarMass	Molar mass of species [g/mol]
// @param initCon		Initial concentration [kg/m^3]
// @param diffCoeff	Diffusion coefficient [m^2/s]
// @param name			Species name
// @param transprot	bool to set if the species is to be tranpsorted with the 
//							fluid velocity
//*****************************************************************************
void meshCell::addSpecies(double molarMass, double initCon, double diffCoeff,
	std::string name, bool transport){

	species spec(molarMass, initCon, diffCoeff, name, transport);
	// loop over connections to see if the species needs to be added
	// to a surface
	for (int conCount = 0; conCount < connections.size(); conCount ++){
		connection* thisCon = getConnection(conCount);
		surface* conSurface = thisCon->getSurface();
		// if the pointer is not null then add species to that surface
		if(conSurface->isInit){
			conSurface->addSpecies(molarMass, initCon, diffCoeff, name, transport);
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
// Gets species concentration [kg/m^3]
//*****************************************************************************
double meshCell::getSpecCon(int specID){
	assert(specID <= speciesVector.size() and specID>= 0);
	species* spec = getSpecies(specID);
	return spec->c;
}

//*****************************************************************************
// Sets the species concentration
//
// @param con		Concentration [kg/m^3]
// @param specID	ID of the species
//*****************************************************************************
void meshCell::setSpeciesConcentration(double con, int specID){
	species spec = speciesVector[specID];
	spec.c = con;
}

//*****************************************************************************
// Sets the cells pressure
//
// @param pressure	Pressure in Pa
//*****************************************************************************
void meshCell::setPressure(double pressure){
	scalarVariables.setPressure(pressure);
}

//*****************************************************************************
// Sets the cells temperature
//
// @param temp		Temperature in kelvin
//*****************************************************************************
void meshCell::setTemperature(double temp){
	scalarVariables.setTemperature(temp);
}

//*****************************************************************************
// Sets the cells scalar neutron flux
//
// @param phi_		Neutron flux in 1/cm^2/s
//*****************************************************************************
void meshCell::setNeutronFlux(double phi_){
	scalarVariables.setNeutronFlux(phi_);
}

//*****************************************************************************
// Sets the cells interfacial area concentration for bubble transport
//
// @param intAreaCon_	Interfacial area concentration 1/m
//*****************************************************************************
void meshCell::setGasInterfacialAreaCon(double intAreaCon_){
	scalarVariables.setGasInterfacialAreaCon(intAreaCon_);
}

//*****************************************************************************
// Sets the cells interfacial area concentration for wall surface area
//
// @param intAreaCon_	Interfacial area concentration 1/m
//*****************************************************************************
void meshCell::setWallInterfacialAreaCon(double intAreaCon_){
	scalarVariables.setWallInterfacialAreaCon(intAreaCon_);
}

//*****************************************************************************
// Sets the cells gas void fraction
//
// @param fract	Gas void fraction
//*****************************************************************************
void meshCell::setGasVoidFraction(double fract){
	scalarVariables.setGasVoidFraction(fract);
}

//*****************************************************************************
// Gets the cells pressure in Pa
//
//*****************************************************************************
double meshCell::getPressure(){
	return scalarVariables.getPressure();
}

//*****************************************************************************
// Gets the cells temperature in kelvin
//
//*****************************************************************************
double meshCell::getTemperature(){
	return scalarVariables.getTemperature();
}

//*****************************************************************************
// Gets the cells scalar neutron flux in 1/cm^2/s
//
//*****************************************************************************
double meshCell::getNeutronFlux(){
	return scalarVariables.getNeutronFlux();
}

//*****************************************************************************
// Gets the cells gas interfacial area concentration in 1/m
//
//*****************************************************************************
double meshCell::getGasInterfacialAreaCon(){
	return scalarVariables.getGasInterfacialAreaCon();
}

//*****************************************************************************
// Gets the cells wall interfacial area concentration in 1/m
//
//*****************************************************************************
double meshCell::getWallInterfacialAreaCon(){
	return scalarVariables.getWallInterfacialAreaCon();
}

//*****************************************************************************
// Gets the cells gas void fraction
//
//*****************************************************************************
double meshCell::getGasVoidFraction(){
	return scalarVariables.getGasVoidFraction();
}

//*****************************************************************************
// Gets a pointer to the scalar data object
//
//*****************************************************************************
scalarData* meshCell::getScalarData(){
	return &scalarVariables;	
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
	for (auto spec : speciesVector){
		spec.clean();
	}
	speciesVector.clear();
	// loop over cell connections
	for (int conCount = 0; conCount < connections.size(); conCount ++){
		connection* thisCon = getConnection(conCount);
	}
}
