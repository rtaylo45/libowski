//*****************************************************************************
// Author: Zack Taylor

// The mesh class. Houses info on the problem domain
//*****************************************************************************
#include <iostream>
#include "modelMesh.h"

//*****************************************************************************
// Builds the problem mesh
//*****************************************************************************
void modelMesh::buildGeometry(){
	
	// Creates the mesh cells
	createCells();
	
	// Creates the cell faces
	createCellFaces();

	// Connects the mesh cells
	connectCells();
}

//*****************************************************************************
// Creates the nodes
//*****************************************************************************
void modelMesh::createCells(){
	int absIndex = 0;

	for (int i = 0; i < numOfxCells; i++){
		for (int j = 0; j < numOfyCells; j++){
			double x = i*dx + dx/2.; // x-coordinate
			double y = j*dy + dy/2.; // y-coordinate

			meshCell cell(i, j, x, y, absIndex);
			meshCells.push_back(cell);	
			absIndex++;
		}
	}
}
//*****************************************************************************
// Creates the cell faces
//*****************************************************************************
void modelMesh::createCellFaces(){
	int absIndex = 0, jmax;
	
	//int k = (numOfxCells+1)*(2*numOfyCells+1) + numOfxCells;
	for (int i = 0; i <= 2*numOfxCells; i++){
		if (i%2){
			jmax = numOfyCells+1;
		}
		else{
			jmax = numOfyCells;
		}
		for (int j = 0; j < jmax; j++){
			meshCellFace face(i, j, absIndex);
			meshCellFaces.push_back(face);
			absIndex++;
		}
	}

}

//*****************************************************************************
// Connects the nodes
//*****************************************************************************
void modelMesh::connectCells(){
	for (int i = 0; i < numOfxCells; i++){
		for (int j = 0; j < numOfyCells; j++){
			meshCell* cell = getCellByLoc(i,j);

			// Sets pointer to the cells
			cell->northCellPtr = getCellByLoc(i,j+1);
			cell->southCellPtr = getCellByLoc(i,j-1);
			cell->westCellPtr = getCellByLoc(i-1,j);
			cell->eastCellPtr = getCellByLoc(i+1,j);

			// Gets the absolute indices for the cell faces
			int kEast = i*(2*numOfyCells+1) + j;
			int kWest = (i+1)*(2*numOfyCells+1) + j;
			int kSouth = numOfyCells*(2*i+1) + i + j;
			int kNorth = numOfyCells*(2*i+1) + i + j + 1;
	
			// Sets the pointers for cell faces	
			cell->eastFacePtr = &meshCellFaces[kEast];
			cell->westFacePtr = &meshCellFaces[kWest];
			cell->southFacePtr = &meshCellFaces[kSouth];
			cell->northFacePtr = &meshCellFaces[kNorth];

			//if (Node->northPtr != nullptr)
			//	std::cout << "North Node " <<Node->northPtr->i << ' ' << Node->northPtr->j << std::endl;
			//if (Node->southPtr != nullptr)
			//	std::cout << "South Node " <<Node->southPtr->i << ' ' << Node->southPtr->j << std::endl;
			//if (Node->westPtr != nullptr)
			//	std::cout << "west Node " <<Node->westPtr->i << ' ' << Node->westPtr->j << std::endl;
			//if (Node->eastPtr != nullptr)
			//	std::cout << "East Node " <<Node->eastPtr->i << ' ' << Node->eastPtr->j << std::endl;
			//std::cout << " " << std::endl;
		}
	}
}

//*****************************************************************************
// Returns a pointer to the node by i,j location
//
// @param i	x index of cell
// @param j y index of cell
//*****************************************************************************
meshCell* modelMesh::getCellByLoc(int i, int j){
	meshCell *returnPtr;

	if (checkCellLoc(i, j)){
		int k = j + i*numOfyCells;
		assert(k <= meshCells.size() and k >= 0);
		returnPtr = &meshCells[k];
	}
	else{
		returnPtr = nullptr;
	}
	return returnPtr;
}

//*****************************************************************************
// Returns a pointer to the node by absolution index
//
// @param k	Absolute index of cell
//*****************************************************************************
meshCell* modelMesh::getCellByLoc(int k){
	meshCell *returnPtr;

	if (k <= meshCells.size() and k >= 0){
		returnPtr = &meshCells[k];
	}
	else{
		returnPtr = nullptr;
	}
	return returnPtr;
}

//*****************************************************************************
// Checks to see if i,j location is valid
//
// @param i	x index of cell
// @param j	y index of cell
//*****************************************************************************
bool modelMesh::checkCellLoc(int i, int j){

	if (i < 0)
		return false;
	else if (i >= numOfxCells)
		return false;
	else if (j < 0)
		return false;
	else if (j >= numOfyCells)
		return false;
	else
		return true;
}

//*****************************************************************************
// Sets a constant x velocity across the whole problem
//
// @param velocity	x velocity [ft/s]
//*****************************************************************************
void modelMesh::setConstantXVelocity(double velocity){
	for (int i = 0; i < numOfxCells; i++){
		for (int j = 0; j < numOfyCells; j++){
			meshCell* cell = getCellByLoc(i,j);
			
			cell->northFacePtr->xVl = velocity;
			cell->southFacePtr->xVl = velocity;
		}
	}
}

//*****************************************************************************
// Sets a constant x velocity across a channel
//
// @param velocity	x velocity [ft/s]
//*****************************************************************************
void modelMesh::setConstantXVelocity(double velocity, int column){
	for (int j = 0; j < numOfyCells; j++){
		meshCell* cell = getCellByLoc(column,j);
		
		cell->northFacePtr->xVl = velocity;
		cell->southFacePtr->xVl = velocity;
	}
}

//*****************************************************************************
// Sets a constant y velocity across the whole problem
//
// @param velocity	y velocity [ft/s]
//*****************************************************************************
void modelMesh::setConstantYVelocity(double velocity){
	for (int i = 0; i < numOfxCells; i++){
		for (int j = 0; j < numOfyCells; j++){
			meshCell* cell = getCellByLoc(i,j);
			
			cell->westFacePtr->yVl = velocity;
			cell->eastFacePtr->yVl = velocity;
		}
	}
}

//*****************************************************************************
// Sets a constant y velocity across a channel
//
// @param velocity	y velocity [ft/s]
//*****************************************************************************
void modelMesh::setConstantYVelocity(double velocity, int row){
	for (int i = 0; i < numOfxCells; i++){
		meshCell* cell = getCellByLoc(i,row);
		
		cell->westFacePtr->yVl = velocity;
		cell->eastFacePtr->yVl = velocity;
	}
}
//*****************************************************************************
// Adds a species to the model
//
// @param molarMass	Molar mass of species [lbm/mol]
// @param [initCon]	Initial concentration [lbm/ft^3]
//*****************************************************************************
int modelMesh::addSpecies(double molarMass, double initCon = 0.0){
	for (int i = 0; i < numOfxCells; i++){
		for (int j = 0; j < numOfyCells; j++){
			meshCell* cell = getCellByLoc(i,j);

			cell->addSpecies(molarMass, initCon);
		}
	}
	int specID = numOfSpecs;
	numOfSpecs++;
	return specID;
}
//*****************************************************************************
// Returns a pointer to the spcies object in a cell
//
//	@param i			x index
// @param j			y index
// @param specID	Species ID
//*****************************************************************************
species* modelMesh::getSpeciesPtr(int i, int j, int specID){
	meshCell* cell = getCellByLoc(i,j);
	species* specPtr = cell->getSpecies(specID);
	return specPtr;
}
//*****************************************************************************
// Gets the species concentration
//
//	@param i			x index
// @param j			y index
// @param specID	Species ID
//*****************************************************************************
double modelMesh::getSpecies(int i, int j, int specID){
	species* spec = getSpeciesPtr(i, j, specID);
	double specCon = spec->c;
	return specCon;
}
//*****************************************************************************
// Sets the source terms for a species in a cell
//
//	@param i			x index
// @param j			y index
// @param specID	Species ID
// @param coeffs	A vector of species coefficients size of number of species
//						[lbm/s]
// @param s			Constant source in cell [lbm/ft^3/s]
//*****************************************************************************
void modelMesh::setSpeciesSource(int i, int j, int specID, std::vector<double> 
		coeffs, double s){
	assert(coeffs.size() == numOfSpecs);
	species* spec = getSpeciesPtr(i, j, specID);
	spec->coeffs = coeffs;
	spec->s = s;
}

//*****************************************************************************
// Cleans the model
//*****************************************************************************
void modelMesh::clean(){

	xLength = 0.0;
	yLength = 0.0;
	numOfxCells = 0;
	numOfyCells = 0;
	numOfTotalCells = 0;
	meshCells.clear();
	meshCellFaces.clear();
	dx = 0.0;
	dy = 0.0;
	cleanSpecies();
}
//*****************************************************************************
// Cleans species in the model
//*****************************************************************************
void modelMesh::cleanSpecies(){
	for (int i = 0; i < numOfxCells; i++){
		for (int j = 0; j < numOfyCells; j++){
			meshCell* cell = getCellByLoc(i,j);
			cell->cleanSpecies();
		}
	}
}
