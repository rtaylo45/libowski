//*****************************************************************************
// Author: Zack Taylor

// The mesh class. Houses info on the problem domain
//*****************************************************************************
#include <iostream>
#include "ModelMeshType.h"

//*****************************************************************************
// Builds the problem mesh
//*****************************************************************************
void modelMesh::buildGeometry(){
	
	// Creates the mesh cells
	createCells();
	
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
	int absIndex = 0;
	
	for (int i = 0; i < numOfxCells+1; i++){
		for (int j = 0; j < numOfyCells+1; j++){
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
			
			cell->northCellPtr = getCellByLoc(i,j+1);
			cell->southCellPtr = getCellByLoc(i,j-1);
			cell->westCellPtr = getCellByLoc(i-1,j);
			cell->eastCellPtr = getCellByLoc(i+1,j);
			//std::cout << "Current Node " << Node->i << ' ' << Node->j << std::endl;
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
	else if (j >= numOfxCells)
		return false;
	else
		return true;
}
