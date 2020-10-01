//*****************************************************************************
// Author: Zack Taylor 
//
//*****************************************************************************
#include "wallDeposition.h"

//**************************************************************************
// Constructor
//
//**************************************************************************
wallDeposition::wallDeposition(){
}

//**************************************************************************
// Sets the coefficients for the model
//
// @param massTransferCoefficient	Mass Transfer coefficient to the wall
//												[m/s]
// @param lID								ID for the liquid species
// @param sID								ID for the surface spcies
// @param infiniteSinkLogic			Logicical for the infinite sink 
//												assumption
//**************************************************************************
void wallDeposition::setModel(double massTransferCoefficient, int lID,
	int sID, bool infiniteSinkLogic){
	massTransferCoeff = massTransferCoefficient;
	liquidID = lID;
	surfaceID = sID;
	infiniteSink = infiniteSinkLogic;
}

//*****************************************************************************
// Gets the coefficient for a species relation
//
// @param otherSpecID		specID of the mass transfer coefficient to find.
// @param scalarVariables	Pointer to scalar data object for the mesh cell
//*****************************************************************************
double wallDeposition::getMassTransferCoeff(int otherSpecID, scalarData* 
	scalarVariables){
	double coeff = 0.0;
	double wallSurfaceAreaCon = scalarVariables->getWallInterfacialAreaCon();

	if (otherSpecID == surfaceID){
		coeff = massTransferCoeff*wallSurfaceAreaCon;
	}
	else if (otherSpecID == liquidID and infiniteSink){
		coeff = -1.*massTransferCoeff*wallSurfaceAreaCon;
	}
	else {
		coeff = 0.0;
	}
	return coeff;
}
