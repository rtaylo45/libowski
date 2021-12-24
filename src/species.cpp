#include <iostream>
#include "species.h"

//**************************************************************************
// Constructor
//
// @param molarMass  Molar mass of species [g/mol]
// @param initCon    Initial concentration [kg/m^3]
// @param diffCoeff  Diffusion coefficient [m^2/s]
// @param name      String of the species name
// @param transport  Bool to set if the species is to be transported by the
//              velocity field
//**************************************************************************
species::species(double molarMass, double initCon, double diffCoeff,
  std::string name_, bool transport_){
  MM = molarMass;
  c = initCon;
  D = diffCoeff;
  name = name_;
  transport = transport_;
}

//**************************************************************************
// Gets the mass transfer coefficient for a species pair
//
// @param otherSpecID  Species ID of the pair combindation
// @param physicsID    Counter of the physcis model
//**************************************************************************
double species::getTransitionCoeff(int otherSpecID, int physicsID,
  scalarData* scalarVariables){
  physicsModel* model = sourceTerms[physicsID];
  return model->getTransitionCoeff(otherSpecID, scalarVariables);
}

//**************************************************************************
// Adds a generic source term to the species
//
// @param coeffs  Vector of coefficients [1/s]
//**************************************************************************
void species::addGenericSourceTerm(std::vector<double> coeffs){
  // Generate model
  physicsModel* model = physicsModelFactory::getPhysicsModel("generic");
  // Set model coefficients
  model->setModel(coeffs);
  // Add model
  addSourceTerm(model);
}

//**************************************************************************
// Adds a neutron induced source term to the species
//
// @param coeffs  Vector of coefficients [cm^2]
//**************************************************************************
void species::addNIRSourceTerm(std::vector<double> coeffs){
  // Generate model
  physicsModel* model = physicsModelFactory::getPhysicsModel("neutronInduced");
  // Set model coefficients
  model->setModel(coeffs);
  // Add model
  addSourceTerm(model);
}

//**************************************************************************
// Adds a wall deposition source term to the species
//
// @param h        MassTransfer coefficient [m/s]
// @param mID      My species ID
// @param lID      ID of the liquid phase species
// @param sID      ID of the surface species
// @param infSink    Logical to set if infinite sink assumption
//**************************************************************************
void species::addWallDepositionSourceTerm(double h, int mID, int lID,
    int sID, bool infSink){
  // Generate model
  physicsModel* model = physicsModelFactory::getPhysicsModel("wallDeposition");
  // Set model coefficients
  model->setModel(h, mID, lID, sID, infSink);
  // Add model
  addSourceTerm(model);
}

//**************************************************************************
// Adds a gas sparging source term to the species
//
// @param k        MassTransfer coefficient [m/s]
// @param H        Henrys law constant [mol/m^3/Pa]
// @param mID      My species ID
// @param lID      ID of the liquid phase species
// @param gID      ID of the surface species
// @param infSink    Logical to set if infinite sink assumption
//**************************************************************************
void species::addGasSpargingSourceTerm(double k, double H, int mID, int lID,
    int gID){
  // Generate model
  physicsModel* model = physicsModelFactory::getPhysicsModel("gasSparging");
  // Set model coefficients
  model->setModel(k, H, MM, mID, lID, gID);
  // Add model
  addSourceTerm(model);
}
//**************************************************************************
// Adds a generic removal term to the species
//
// @param coeff    Removal rate [1/s]
// @param mID      My species ID
//**************************************************************************
void species::addGenericRemovalSourceTerm(double coeff, int mID){
  // Generate model
  physicsModel* model = physicsModelFactory::getPhysicsModel("genericRemoval");
  // Set model coefficients
  model->setModel(mID, coeff);
  // Add model
  addSourceTerm(model);
}

//**************************************************************************
// Addes a source term model
//
// @param sourceTermModel  Pointer to physics source term model
//**************************************************************************
void species::addSourceTerm(physicsModel* sourceTermModel){
  sourceTerms.push_back(sourceTermModel);
}

//**************************************************************************
// Clean
//**************************************************************************
void species::clean(){
  c = 0.0;
  MM = 0.0;
  s = 0.0;
  D = 0.0;
  name = "None";
  for (auto p : sourceTerms){
    delete p;
  }
  sourceTerms.clear();
}
