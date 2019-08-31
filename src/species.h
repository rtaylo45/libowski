//*****************************************************************************
// Author: Zack Taylor
//
// Defines a chemical or iostopic conserved species
//*****************************************************************************
#include <vector>

class species {

	// class attributes
	public:
	// Concentration [lbm/ft^3]
	double c = 0.0;
	// molar mass [lbm/mol]
	double MM = 0.0;
	// Constant volumetric source terms [lbm/ft^3/s]
	double s = 0.0;
	// Vector of linear source term coefficients in order of species IDs
	std::vector<double> coeffs;

	// Class methods
	public:
	// Constructor
	species(double, double);

	// Clean
	void clean();
};
