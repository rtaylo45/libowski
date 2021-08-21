//*****************************************************************************
// Author: Zack Taylor
//
// Cell face type. This class defines the sides of each mesh data cell. The
// flux of species going into or out of each mesh cell is defined by the 
// velocity of the fluid accross the cell faces. Each cell face will have two
// owners. If the face is parallel to the x-axis then the owners will be a
// mesh cell below and a mesh cell above. If the face is parallel to the y-axis
// then the mesh cells will be on the righ and left. Positive flow in the y
// direction is defined to go from bottom to top. Positive flow in the x
// direction is defined to go from left to right.
//*****************************************************************************
#ifndef MESHCELLFACE_H
#define MESHCELLFACE_H
#include "surface.h"

class meshCellFace {

	public:
	// Index of face in x direction
	int i;
	// Index of face in y direction
	int j;
	// Absolute index of the surface
	int absIndex;
	// Velocity component normal to the surface
	double vl = 0.0;
	// x direction length [m]
	double dx = 0.0;
	// y direction length [m]
	double dy = 0.0;
	// Surface object
	surface mySurface = surface();

	// Constructor
	meshCellFace(int, int, int);
	// Gets a pointer to the surface
	surface* getSurface();

};
#endif
