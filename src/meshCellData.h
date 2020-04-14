//*****************************************************************************
// Author: Zack Taylor
//
// Mesh cell data type. This class holds scalar cell data information based
// on a 2-D finite volume discritization 
//*****************************************************************************
#ifndef MESHCELLDATA_H
#define MESHCELLDATA_H
#include "meshCellFace.h"
#include "species.h"
#include <vector>
#include <assert.h>
#include <string>

class meshCell {

	// Class attributes
	public:
	// Cell index in the x direction
	int i = -1;
	// Cell index in the y direction
	int j = -1;
	// Absolut index
	int absIndex = -1;
	// X position defined at the center of the cell
	double x = -1.;
	// Y position defined at the cent of the cell
	double y = -1. ;
	// dx of cell
	double dx = -1.;
	// dy of cell
	double dy = -1.;
	// Temperature of cell in kelvin
	double T = -1.;
	// Presure in lbf/in^2
	double P = -1.;
	// Denote if a cell is at a  boundary
	bool boundary = false;
	// Boundary location: 0 = North
	//							 1 = South
	//							 2 = East
	//							 3 = West
	int boundaryLoc = -1;
	// Boundary condition type
	std::string boundaryType = "None";
	// Vector of the species in the cell
	std::vector<species> speciesVector;

	// Connection info for mesh cells
	// Pointer to the east mesh cell
	meshCell *eastCellPtr = nullptr;
	// Pointer to the west mesh cell
	meshCell *westCellPtr = nullptr;
	// Pointer to the north mesh cell
	meshCell *northCellPtr = nullptr;
	// Pointer to the south mesh cell
	meshCell *southCellPtr = nullptr;

	// Connection info for mesh cell faces
	// Pointer to the east mesh cell
	meshCellFace *eastFacePtr = nullptr;
	// Pointer to the west mesh cell
	meshCellFace *westFacePtr = nullptr;
	// Pointer to the north mesh cell
	meshCellFace *northFacePtr = nullptr;
	// Pointer to the south mesh cell
	meshCellFace *southFacePtr = nullptr;

	public:
	// Constructor
	meshCell(int, int, int, double, double, double, double);
	// Add species
	void addSpecies(double, double, double);
	// Gets a pointer to the species
	species* getSpecies(int);
	// Gets this cells species concentration
	double getSpecCon(int);
	// Sets species concentration
	void setSpeciesConcentration(double, int);
	// Clean species
	void cleanSpecies();


};
#endif
