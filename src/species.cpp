#include <iostream>
#include "species.h"

//**************************************************************************
// Constructor
//
// @param molarMass	Molar mass of species [g/mol]
// @param initCon		Initial concentration [kg/m^3]
// @param diffCoeff	Diffusion coefficient [m^2/s]
// @param name			String of the species name
// @param transport	Bool to set if the species is to be transported by the 
//							velocity field
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
// Addes a row of source term coefficients to the source term array
//
// @param coeffs	Array of coefficients 
//**************************************************************************
void species::addCoeffRow(std::vector<double> v){
	ArrayD coeffRow;
	coeffRow = Eigen::Map<Eigen::ArrayXd>(v.data(), v.size()).transpose();	
	if (coeffs.rows() != 0 and coeffs.cols() != 0){
		assert(coeffs.cols() == coeffRow.cols() and coeffRow.rows() == 1);
	}
	// Coefficient array has not been set yet
	if (coeffs.rows() == 0 and coeffs.cols() == 0){
		coeffs = coeffRow;
	}
	// Coefficient array already has a row of values
	else{
		// Define the new coefficient array with 1+ row
		ArrayD newCoeffs = ArrayD::Zero(coeffs.rows()+1, coeffs.cols());
		// Sets the new coefficient array with the previous coefficients
		newCoeffs.topRows(coeffs.rows()) = coeffs;
		// Adds the new coeff row to the last row of the new coeff array
		newCoeffs.row(coeffs.rows()) = coeffRow;
		// Sets the coefficient array to the new coefficients
		coeffs = newCoeffs;
	}
}
//**************************************************************************
// Gets the mass transfer coefficient for a species pair
//
// @param otherSpecID	Species ID of the pair combindation
// @param physicsID		Counter of the physcis model
//**************************************************************************
double species::getMassTransferCoeff(int otherSpecID, int physicsID, 
	scalarData* scalarVariables){
	physicsModel* model = sourceTerms[physicsID];
	return model->getMassTransferCoeff(otherSpecID, scalarVariables);
}

//**************************************************************************
// Adds a generic source term to the species
//
// @param coeffs	Vector of coefficients [1/s]
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
// @param coeffs	Vector of coefficients [cm^2]
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
// Addes a source term model 
//
// @param sourceTermModel	Pointer to physics source term model
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
	coeffs = ArrayD();
	for (auto p : sourceTerms){
		delete p;
	}
	sourceTerms.clear();
}
