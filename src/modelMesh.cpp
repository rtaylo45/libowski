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
			//std::cout << i << " " << j << std::endl;
			meshCell* cell = getCellByLoc(i,j);
			//std::cout << cell<< std::endl;
			//std::cout << "Current Node " << cell->i << ' ' << cell->j << std::endl;
			cell->northCellPtr = getCellByLoc(i,j+1);
			cell->southCellPtr = getCellByLoc(i,j-1);
			cell->westCellPtr = getCellByLoc(i-1,j);
			cell->eastCellPtr = getCellByLoc(i+1,j);

			// Hope i got this part right
			int kEast = i*(2*numOfyCells+1) + j;
			int kWest = (i+1)*(2*numOfyCells+1) + j;
			int kSouth = numOfyCells*(2*i+1) + i + j;
			int kNorth = numOfyCells*(2*i+1) + i + j + 1;
		
			//std::cout << i << " " << j << " " << kEast << " " << kWest << " " <<
			//	kSouth << " " << kNorth << std::endl;	
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
//*****************************************************************************
void modelMesh::setConstantYVelocity(double velocity){
	for (int i = 0; i < numOfxCells; i++){
		for (int j = 0; j < numOfyCells; j++){
			meshCell* cell = getCellByLoc(i,j);
			
			cell->northFacePtr->yVl = velocity;
			cell->southFacePtr->yVl = velocity;
		}
	}
}

//*****************************************************************************
// Sets a constant x velocity across a channel
//*****************************************************************************
void modelMesh::setConstantYVelocity(double velocity, int row){
	for (int i = 0; i < numOfxCells; i++){
		meshCell* cell = getCellByLoc(i,row);
		
		cell->northFacePtr->yVl = velocity;
		cell->southFacePtr->yVl = velocity;
	}
}

