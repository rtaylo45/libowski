#include "cellConnection.h"

//**************************************************************************
// Constructor
//
// @param conCell			Cell at other end of connection
// @param conFace			The surface between the two connecting cells
// @param area_			The area of the surface that the fluid is flowing 
//								through. This will be orthogonal to the velocity
//								vector which is normal to the surface. Because
//								libowski is 2D, the area is just dx or dy.
// @param direction_		Direction indecator for the convection transition
//								coefficient
//								nTran = -1
//								sTran = 1
//								eTran = -1
//								wTran = 1
//	@param loc				Index for the location of the connection face
//								north = 0
//								south = 1
//								east = 2
//								west = 3
//	@param distance		The distance from the principal cell center to the 
//								connection cell center. These are the small delta
//								y and x from my PhD thesis.
//**************************************************************************
connection::connection(meshCell *conCell, meshCellFace *conFace, double 
	area_, double direction_, int loc_, double distance_){

	// Set connection variables
	connectionCellPtr = conCell;
	connectionFacePtr = conFace;
	area = area_;
	direction = direction_;
	loc = loc_;
	distance = distance_;
}
