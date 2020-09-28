//*****************************************************************************
// Author: Zack Taylor
//
// Mesh cell data type. This class holds scalar cell data information based
// on a 2-D finite volume discritization 
//*****************************************************************************
#ifndef MESHCELLDATA_H
#define MESHCELLDATA_H
#include <vector>
#include <iostream>
#include <assert.h>
#include <string>

#include "cellConnection.h"
#include "meshCellFace.h"
#include "surface.h"
#include "species.h"
#include "scalarData.h"

// Forward decleration
class connection;

class meshCell {

	// Class attributes
	public:
	// Cell index in the x direction
	int i = -1;
	// Cell index in the y direction
	int j = -1;
	// Absolut index
	int absIndex = -1;
	// X position defined at the center of the cell [m]
	double x = -1.;
	// Y position defined at the cent of the cell [m]
	double y = -1. ;
	// dx of cell [m]
	double dx = -1.;
	// dy of cell [m]
	double dy = -1.;
	// cell volume [m^2] because its 2D
	double volume = 0.0;
	// Vector of cell connections
	std::vector<connection> connections;
	// Flag to set if the second order convective upwind flux is used
	bool secondOrderFlux = true;

	private:
	// Vector of the species in the cell
	std::vector<species> speciesVector;
	// Scalar data object
	scalarData scalarVariables;

	public:
	// Constructor
	meshCell(int, int, int, double, double, double, double);
	// Add species
	void addSpecies(double, double = 0.0, double = 0.0, std::string = "None",
		bool = true);
	// Gets a pointer to the species
	species* getSpecies(int);
	// Gets this cells species concentration
	double getSpecCon(int);
	// Sets species concentration
	void setSpeciesConcentration(double, int);
	// Sets the cells temperature
	void setTemperature(double);
	// Sets the cells pressure
	void setPressure(double);
	// Sets the cells scalar neutron flux
	void setNeutronFlux(double);
	// Sets the cells interfacial area concentration 
	void setInterfacialAreaCon(double);
	// Gets the cells temperature
	double getTemperature();
	// Gets the cells pressure
	double getPressure();
	// Gets the cells scalar neutron flux
	double getNeutronFlux();
	// Gets the cells interfacial area concentration 
	double getInterfacialAreaCon();
	// Gets a pointer to the scalar data object
	scalarData* getScalarData();
	// Gets a pointer to the connection
	connection* getConnection(int);
	// Adds a surface 
	void addSurface(int);
	// Clean species
	void cleanSpecies();

};
#endif
