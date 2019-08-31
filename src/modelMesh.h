//*****************************************************************************
// Author: Zack Taylor
// 
// The model mesh. Houses info on the problem domain
//*****************************************************************************
#ifndef MODELMESH_H
#define MODELMESH_H
#include "meshCellData.h"
#include "meshCellFace.h"
#include <vector>
#include <assert.h>

class modelMesh {

	public:
	// Total length in x direction
	double xLength = 0.0;
	// Total length in y direction
	double yLength = 0.0;
	// Number of cells in the x direction
	int numOfxCells = 0;
	// Number of cells in the y direction
	int numOfyCells = 0;
	// Total number of cells
	int numOfTotalCells = 0;
	// Vector of all cells
	std::vector<meshCell> meshCells;
	// Vector of all cell faces
	std::vector<meshCellFace> meshCellFaces;
	// Change in x direction
	double dx = 0.0;
	// Change in y direction
	double dy = 0.0;
	// Number of specie in the model
	int numOfSpecs = 0;
	
	public:
	// Constructor
	modelMesh(int, int, double, double);
	// Builds the geometry
	void buildGeometry();
	// Gets the node by location from i,j
	meshCell* getCellByLoc(int, int);
	// Gets the node by location from absolution index k
	meshCell* getCellByLoc(int);
	// Sets a constant x velocity across the whole problem
	void setConstantXVelocity(double);
	// Sets a constant x velocity across a column of cells
	void setConstantXVelocity(double, int);
	// Sets a constant y velocity across the whole problem
	void setConstantYVelocity(double);
	// Sets a constant y velocity across a row of cells
	void setConstantYVelocity(double, int);
	// Adds a species to the system
	int addSpecies(double, double);
	// Gets a pointer to the spcies object
	species* getSpeciesPtr(int, int, int);
	// Gets the species concentration
	double getSpecies(int, int, int);
	// Sets the species source terms
	void setSpeciesSource(int, int, int, std::vector<double>, double);
	// Cleans the model
	void clean();
	// Cleans species
	void cleanSpecies();

	private:
	// Creates the cells
	void createCells();
	// Connects the cells
	void connectCells();
	// Creates the cell faces
	void createCellFaces();
	// Checks the i,j for validity
	bool checkCellLoc(int, int);

};
#endif
