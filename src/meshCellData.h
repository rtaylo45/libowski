//*****************************************************************************
// Author: Zack Taylor
//
// Mesh cell data type. This class holds scalar cell data information based
// on a 2-D finite volume discritization 
//*****************************************************************************
#include "meshCellFace.h"
#include "species.h"
#include <vector>

class meshCell {

	// Class attributes
	public:
	// Cell index in the x direction
	int i;
	// Cell index in the y direction
	int j;
	// Absolut index
	int absIndex;
	// X position defined at the center of the cell
	double x;
	// Y position defined at the cent of the cell
	double y;
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
	//**************************************************************************
	// Constructor
	//
	// @param iIndex			Index of cell in x direction
	// @param jIndex			Index of cell in y direction
	// @param absoluteIndex	Absolute index of the cell
	// @param xCor				Location of cell center in x direction
	// @param yCor				Location of cell center in y direction
	//**************************************************************************
	meshCell(int iIndex, int jIndex, int absoluteIndex, double xCor, double yCor){
		i = iIndex;
		j = jIndex;
		absIndex = absoluteIndex;
		x = xCor;
		y = yCor;	
	}
	//**************************************************************************
	// Add species
	//**************************************************************************
	void addSpecies(double, double);

};
