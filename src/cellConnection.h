//*****************************************************************************
// Author: Zack Taylor
//
// Connection class between cells
//*****************************************************************************
#ifndef CELLCONNECTION_H
#define CELLCONNECTION_H
#include "meshCellData.h"
#include "meshCellFace.h"

// Forward decleration
class meshCell;
class meshCellFace;

class connection {

	public:
	// Cross sectional area for transport through the connection
	double area = 0.0;
	// direction used to define the default flow direction for libowski
	double direction = 0.0;
	// The cell face that is between the two connecting cells
	meshCellFace* connectionFacePtr = nullptr;
	// pointer to connecting cell
	meshCell* connectionCellPtr = nullptr;

	// constructor
	connection(meshCell*, meshCellFace*, double, double);
};

#endif
