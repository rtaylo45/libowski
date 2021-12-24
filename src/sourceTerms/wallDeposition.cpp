//*****************************************************************************
// Author: Zack Taylor
//
//*****************************************************************************
#include "wallDeposition.h"
#include <iostream>

//**************************************************************************
// Constructor
//
//**************************************************************************
wallDeposition::wallDeposition(){
}

//**************************************************************************
// Sets the coefficients for the model
//
// @param massTransferCoefficient  Mass Transfer coefficient to the wall
//                        [m/s]
//  @param mID                ID for spec that owns me
// @param lID                ID for the liquid species
// @param sID                ID for the surface spcies
// @param infiniteSinkLogic      Logicical for the infinite sink
//                        assumption
//**************************************************************************
void wallDeposition::setModel(double massTransferCoefficient, int mID,
  int lID, int sID, bool infiniteSinkLogic){
  massTransferCoeff = massTransferCoefficient;
  myID = mID;
  liquidID = lID;
  surfaceID = sID;
  infiniteSink = infiniteSinkLogic;
}

//*****************************************************************************
// Gets the coefficient for a species relation
//
// @param otherSpecID    specID of the mass transfer coefficient to find.
// @param scalarVariables  Pointer to scalar data object for the mesh cell
//*****************************************************************************
double wallDeposition::getTransitionCoeff(int otherSpecID, scalarData*
  scalarVariables){
  double coeff = 0.0;
  double wallSurfaceAreaCon = scalarVariables->getWallInterfacialAreaCon();

  //  dC_surface/dt = (h*A/V)*(C_liquid - C_surface)
  if (myID == surfaceID){
    if (otherSpecID == liquidID){
      coeff = massTransferCoeff*wallSurfaceAreaCon;
    }
    else if (otherSpecID == surfaceID and infiniteSink){
      coeff = -1.*massTransferCoeff*wallSurfaceAreaCon;
    }
  }
  else if (myID == liquidID){
    if (otherSpecID == liquidID){
      coeff = -1.*massTransferCoeff*wallSurfaceAreaCon;
    }
    else if (otherSpecID == surfaceID and infiniteSink){
      coeff = massTransferCoeff*wallSurfaceAreaCon;
    }
  }
  return coeff;
}
