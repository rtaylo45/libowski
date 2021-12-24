#include "meshCellFace.h"

//**************************************************************************
// Constructor
//
// @param iIndex      x-direction mesh index
// @param jIndex      y-direction mesh index
// @param absoluteIndex  Total mesh index
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
