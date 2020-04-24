//*****************************************************************************
// Author: Zack Taylor
//
// Connection class between cells
//*****************************************************************************
#ifndef CELLCONNECTION_H
#define CELLCONNECTION_H
#include "meshCellData.h"
#include "meshCellFace.h"
#include <string>

// Forward decleration
class meshCell;
class meshCellFace;

class connection {

	public:
	// Cross sectional area for transport through the connection
	double area = 0.0;
	// direction used to define the default flow direction for libowski
	double direction = 0.0;
	// location of the face connection
	//							 0 = North
	//							 1 = South
	//							 2 = East
	//							 3 = West
	int loc = -1;
	// The distance from the principal cell cent to the connection
	// cell center.
	double distance = 0.0;
	// Denote if a connection is at a  boundary
	bool boundary = false;
	// Boundary condition type
	std::string boundaryType = "None";
	// The cell face that is between the two connecting cells
	meshCellFace* connectionFacePtr = nullptr;
	// pointer to connecting cell
	meshCell* connectionCellPtr = nullptr;
	// Adds a surface to the connection face pointer
	void addSurface();
	// Gets a pointer to the connection surface
	surface* getSurface();

	// constructor
	connection(meshCell*, meshCellFace*, double, double, int, double);
};

#endif
