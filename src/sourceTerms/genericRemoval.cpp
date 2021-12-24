//*****************************************************************************
// Author: Zack Taylor
//
//  Defines the generic source term model
//*****************************************************************************
#include "genericRemoval.h"
#include <iostream>

//**************************************************************************
// Constructor
//
//**************************************************************************
genericRemoval::genericRemoval(){
}

//*****************************************************************************
// Sets the coefficients for the model
//
//  @param myID_  id of the spcies that ones me
// @param coeff  Generic removal source term [1/s]
//*****************************************************************************
void genericRemoval::setModel(int myID_, double coeff){
  myID = myID_;
  removalCoeff = coeff;
}

//*****************************************************************************
// Gets the coefficient for a species relation
//
// @param otherSpecID    specID of the mass transfer coefficient to find.
// @param scalarVariables  Pointer to scalar data object for the mesh cell
//*****************************************************************************
double genericRemoval::getTransitionCoeff(int otherSpecID, scalarData* scalarVariables){
  if (myID == otherSpecID){
    return removalCoeff;
  }

  else {
    return 0.0;
  }
}
