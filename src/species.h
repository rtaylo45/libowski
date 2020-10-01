//*****************************************************************************
// Author: Zack Taylor
//
// Defines a chemical or iostopic conserved species
//*****************************************************************************
#ifndef SPECIES_H
#define SPECIES_H
#include <vector>
#include <string>
#include "matrixTypes.h"
#include "physicsModels.h"
#include "scalarData.h"

class species {

	// class attributes
	public:
	// Name of the species
	std::string name = "None";
	// Concentration [kg/m^3]
	double c = 0.0;
	// molar mass [g/mol]
	double MM = 0.0;
	// Constant volumetric source terms [kg/m^3/s]
	double s = 0.0;
	// Holds the boundary condition value [kg/m^3] or [kg/m^2]
	double bc = 0.0;
	// Diffusion coefficient [m^2/s]
	double D = 0.0;
	// Bool to set if the species is transported with the fluid velocity
	bool transport = true;
	// Vector that holds the pointers to the physics models. These models
	// generate source term coefficients for the transition matrix
	std::vector<physicsModel*> sourceTerms;

	// Class methods
	public:
	// Constructor
	species(double, double = 0.0, double = 0.0, std::string = "None", 
		bool = true);
	// Gets the mass transfer coefficient 
	double getTransitionCoeff(int, int, scalarData*);
	// Adds generic source term
	void addGenericSourceTerm(std::vector<double>);
	// Adds neutron induced reactions source term
	void addNIRSourceTerm(std::vector<double>);
	// Add wall deposition source term
	void addWallDepositionSourceTerm(double, int, int, bool);
	// Clean
	void clean();

	private:
	// Adds a physics model to the source term coefficient vector
	void addSourceTerm(physicsModel*);
	
};
#endif
