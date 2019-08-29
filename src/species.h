//*****************************************************************************
// Author: Zack Taylor
//
// Defines a chemical or iostopic conserved species
//*****************************************************************************
class species {

	// class attributes
	public:
	// Concentration [lbm/ft^3]
	double c = 0.0;
	// molar mass [lbm/mol]
	double MM = 0.0;

	// Class methods
	public:
	void clean(){
		c = 0.0;
		MM = 0.0;
	}

};
