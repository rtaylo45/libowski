//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include "scalarData.h"

//*****************************************************************************
// Constructor
//*****************************************************************************
scalarData::scalarData(){};

//*****************************************************************************
// Sets temperature
//
// @param temp		Temperature in Kelvin
//*****************************************************************************
void scalarData::setTemperature(double temp){
	assert(temp >= 0.0);
	T = temp;
}

//*****************************************************************************
// Sets pressure
//
// @param pres		Pressure in Pa
//*****************************************************************************
void scalarData::setPressure(double pres){
	assert(pres >= 0.0);
	P = pres;
}

//*****************************************************************************
// Sets neutron flux
//
// @param flux		Neutron flux in 1/cm^2/s
//*****************************************************************************
void scalarData::setNeutronFlux(double flux){
	assert(flux >= 0.0);
	phi = flux;
}

//*****************************************************************************
// Sets the interfacial area concentration
//
// @param a		Interfacial area concentration 1/m
//*****************************************************************************
void scalarData::setInterfacialAreaCon(double a){
	assert(a >= 0.0);
	intAreaCon = a;
}

//*****************************************************************************
// Gets temperature in Kelvin
//
//*****************************************************************************
double scalarData::getTemperature(){
	assert(T >= 0.0);
	return T;
}

//*****************************************************************************
// Gets pressure in Pa
//
//*****************************************************************************
double scalarData::getPressure(){
	assert(P >= 0.0);
	return P;
}

//*****************************************************************************
// Gets neutron flux in 1/cm^2/s
//
//*****************************************************************************
double scalarData::getNeutronFlux(){
	assert(phi >= 0.0);
	return phi;
}
//*****************************************************************************
// Gets the interfacial area concentration in 1/m
//
//*****************************************************************************
double scalarData::getInterfacialAreaCon(){
	assert(intAreaCon >= 0.0);
	return intAreaCon;
}
