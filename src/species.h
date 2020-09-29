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
	// Array that holds linear source term coefficients. The number of rows 
	// is equal to the number of species in the system. The number of columns
	// is equal to the number of physics source term models in the system. 
	// Each column represents a physical model, with coefficients for each 
	// species in a row column index. Because difference physics can have
	// different source term models the units in the coefficient array can
	// very, but will be constant along rows. If the physics model is 
	// generic the coefficient row will have units of [1/s]. If the physics
	// is for neutron induced reactions, the coefficients will have units
	// of [cm^2] and thus will be multiplied by the neutron flux before
	// being added to the transition matrix. The neutron flux is housed in
	// the mesh cell data object. These are only two examples. It is done
	// in this way to separate terms in the sourse term models that are 
	// constant and those that very as a function of time or space. 
	ArrayD coeffs;
	// Vector that holds the pointers to the physics models. These models
	// generate source term coefficients for the transition matrix
	std::vector<physicsModel*> sourceTerms;

	// Class methods
	public:
	// Constructor
	species(double, double = 0.0, double = 0.0, std::string = "None", 
		bool = true);
	// Adds a row of coeffs to the source term array
	void addCoeffRow(std::vector<double>);
	// Adds a physics model to the source term coefficient vector
	void addSourceTerm(physicsModel*);

	// Clean
	void clean();
};
#endif
