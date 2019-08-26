//*****************************************************************************
// Author: Zack Taylor
// 
// The model mesh. Houses info on the problem domain
//*****************************************************************************
#include "MeshCellDataType.h"
#include <vector>
#include <assert.h>

class modelMesh {

	public:
	// Total length in x direction
	double xLength;
	// Total length in y direction
	double yLength;
	// Number of cells in the x direction
	int numOfxCells;
	// Number of cells in the y direction
	int numOfyCells;
	// Total number of cells
	int numOfTotalCells;
	// Vector of all cells
	std::vector<meshCell> meshCells;
	// Vector of all cell faces
	std::vector<meshCellFace> meshCellFaces;
	// Change in x direction
	double dx;
	// Change in y direction
	double dy;
	
	public:
	//**************************************************************************
	// Constructor
	// 
	// @param xCells		Number of cells in the x direction
	// @param yCells		Number of cells in the y direction
	// @param xDirLength	Total length in the x direction
	// @param yDirLength	Total length in the y direction
	//**************************************************************************
	modelMesh(int xCells, int yCells, double xDirLength, double yDirLength){
		numOfxCells = xCells;
		numOfyCells = yCells;
		xLength = xDirLength;
		yLength = yDirLength;
		dx = xLength/(float)numOfxCells;
		dy = yLength/(float)numOfyCells;
		numOfTotalCells = numOfxCells*numOfyCells;

		// Builds the geometry
		buildGeometry();
	}

	// Builds the geometry
	void buildGeometry();
	// Gets the node by location from i,j
	meshCell* getCellByLoc(int, int);
	// Gets the node by location from absolution index k
	meshCell* getCellByLoc(int);

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
