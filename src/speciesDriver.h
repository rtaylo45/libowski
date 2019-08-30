//*****************************************************************************
// Author: Zack Taylor 
//
// Driver class for transported species. This class is used to build the 
// species transport problem
//*****************************************************************************
#include "modelMesh.h"

class speciesDriver {

	// Class attributes
	public:
	// Pointer to the model object
	modelMesh* modelPtr = nullptr;
	// Number of species in the model
	int numOfSpecs = 0;

	// Class methods
	public:
	//**************************************************************************
	// Constructor
	//
	// @param model A pointer to the mesh model
	//**************************************************************************
	speciesDriver(modelMesh* model){
		modelPtr = model;
	}
};
