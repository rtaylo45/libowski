//*****************************************************************************
// Author: Zack Taylor
//
// Defines a chemical or iostopic conserved species
//*****************************************************************************
#ifndef SPECIES_H
#define SPECIES_H
#include <vector>
#include <string>

class species {

	// class attributes
	public:
	// Name of the species
	std::string name = "None";
	// Concentration [lbm/ft^3]
	double c = 0.0;
	// molar mass [lbm/mol]
	double MM = 0.0;
	// Constant volumetric source terms [lbm/ft^3/s]
	double s = 0.0;
	// Holds the boundary condition value [lbm/ft^3]
	double bc = 0.0;
	// Diffusion coefficient
	double D = 0.0;
	// Vector of linear source term coefficients in order of species IDs [1/s]
	std::vector<double> coeffs;
	// Vector of linear source terms for neturon induces reactions in order of 
	// species IDs [ft^2]. These values need to be multiplied by the scalar 
	// neutron flux in the cell.
	std::vector<double> transCoeffs;

	// Class methods
	public:
	// Constructor
	species(double, double = 0.0, double = 0.0, std::string = "None");

	// Clean
	void clean();
};
#endif
