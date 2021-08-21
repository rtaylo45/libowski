//*****************************************************************************
// Author: Zack Taylor 
//
//*****************************************************************************
#include "gasSparging.h"

//**************************************************************************
// Constructor
//
//**************************************************************************
gasSparging::gasSparging(){
}

//**************************************************************************
// Sets the coefficients for the model
//
// @param massTransferCoefficient	Mass Transfer coefficient to the gas
//												[m/s]
//	@param Hlc								Henrys law constant [mol/m^3/Pa]
//	@param mmin								Species molar mass [g/mol]
//	@param mID								ID for spec that owns me
// @param lID								ID for the liquid species
// @param gID								ID for the gas species
//**************************************************************************
void gasSparging::setModel(double massTransferCoefficient, double Hlc, double mmin, 
	int mID, int lID, int gID){
	assert(Hlc >= 0.0);
	massTransferCoeff = massTransferCoefficient;
	H = Hlc;
	mm = mmin;
	myID = mID;
	liquidID = lID;
	gasID = gID;
}

//*****************************************************************************
// Gets the coefficient for a species relation
//
// @param otherSpecID		specID of the mass transfer coefficient to find.
// @param scalarVariables	Pointer to scalar data object for the mesh cell
//*****************************************************************************
double gasSparging::getTransitionCoeff(int otherSpecID, scalarData* 
	scalarVariables){
	double coeff = 0.0;
	double gasSurfaceAreaCon = scalarVariables->getGasInterfacialAreaCon();
	double gasVoidFraction = scalarVariables->getGasVoidFraction();
	double liqVoidFraction = 1. - gasVoidFraction;
	double temp = scalarVariables->getTemperature();

	if (myID == gasID){
		// (h*A/V)*1/(liq_void_fraction)
		if (otherSpecID == liquidID){
			coeff = massTransferCoeff*gasSurfaceAreaCon/liqVoidFraction;
		}
		//	-(h*A/V)*1/(gas_void_fraction*molar_mass)*1000*H*R*T
		else if (otherSpecID == gasID){
			coeff = -1.*massTransferCoeff*gasSurfaceAreaCon*1000*H*idealGasR*temp/(gasVoidFraction*mm);
		}
	}
	else if (myID == liquidID){
		//	-(h*A/V)*1/(liq_void_fraction)
		if (otherSpecID == liquidID){
			coeff = -1.*massTransferCoeff*gasSurfaceAreaCon/liqVoidFraction;
		}
		// (h*A/V)*1/(gas_void_fraction*molar_mass)*1000*H*R*T
		else if (otherSpecID == gasID){
			coeff = massTransferCoeff*gasSurfaceAreaCon*1000*H*idealGasR*temp/(gasVoidFraction*mm);
		}
	}
	return coeff;
}
