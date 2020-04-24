#include "meshCellFace.h"

//**************************************************************************
// Constructor
//**************************************************************************
meshCellFace::meshCellFace(int iIndex, int jIndex, int absoluteIndex){
	i = iIndex;	
	j = jIndex;	
	absIndex = absoluteIndex;	
}

//**************************************************************************
// Returns a pointer to the surface object
//**************************************************************************
surface* meshCellFace::getSurface(){
	return &mySurface;
}
