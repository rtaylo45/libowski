#include "cellConnection.h"

//**************************************************************************
// Constructor
//**************************************************************************
connection::connection(meshCell *conCell, meshCellFace *conFace, double 
	area_, double direction){

	// Set connection variables
	connectionCellPtr = conCell;
	connectionFacePtr = conFace;
	area = area_;
	direction = direction_;
}
