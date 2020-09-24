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
	// Vector of linear source term coefficients in order of species IDs [1/s]
	std::vector<double> coeffs;
	// Vector of linear source terms for neturon induces reactions in order of 
	// species IDs [cm^2]. These values need to be multiplied by the scalar 
	// neutron flux in the cell.
	std::vector<double> transCoeffs;
	// Array that holds linear source term coefficients. The number of rows 
	// is equal to the number of species in the system. The number of columns
	// is a minimal of 1. The second row will hold coefficients for neutron 
	// induced reactions [cm^2]. These values need to be multiplied by the 
	// scalar neutron flux before 
	// The frist row holds
	// generice source coefficients which the user can define. 
	ArrayD coeffs;

	// Class methods
	public:
	// Constructor
	species(double, double = 0.0, double = 0.0, std::string = "None", 
		bool = true);

	// Clean
	void clean();
};
#endif
