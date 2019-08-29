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
	//**************************************************************************
	// Constructor
	//
	// @param molarMass	Molar mass of species [lbm/mol]
	// @param initCon		Initial concentration [lbm/ft^3]
	//**************************************************************************
	species(double molarMass, double initCon){
		MM = molarMass;
		c = initCon;
	}

	//**************************************************************************
	// Clean
	//**************************************************************************
	void clean(){
		c = 0.0;
		MM = 0.0;
	}
};
