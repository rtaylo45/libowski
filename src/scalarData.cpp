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
// Sets the interfacial area concentration for gas phase
//
// @param a		Interfacial area concentration 1/m
//*****************************************************************************
void scalarData::setGasInterfacialAreaCon(double a){
	assert(a >= 0.0);
	gasIntAreaCon = a;
}

//*****************************************************************************
// Sets the interfacial area concentration for wall surface
//
// @param a		Interfacial area concentration 1/m
//*****************************************************************************
void scalarData::setWallInterfacialAreaCon(double a){
	assert(a >= 0.0);
	wallIntAreaCon = a;
}

//*****************************************************************************
// Sets the gas void fraction
//
// @param fract		Gas void fraction [fraction]
//*****************************************************************************
void scalarData::setGasVoidFraction(double fract){
	assert(fract >= 0.0);
	gasVoidFraction = fract;
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
// Gets the interfacial area concentration of gas phase in 1/m
//
//*****************************************************************************
double scalarData::getGasInterfacialAreaCon(){
	assert(gasIntAreaCon >= 0.0);
	return gasIntAreaCon;
}
//*****************************************************************************
// Gets the interfacial area concentration wall surface in 1/m
//
//*****************************************************************************
double scalarData::getWallInterfacialAreaCon(){
	assert(wallIntAreaCon >= 0.0);
	return wallIntAreaCon;
}
//*****************************************************************************
// Gets the gas void fraction [fraction]
//
//*****************************************************************************
double scalarData::getGasVoidFraction(){
	assert(gasVoidFraction >= 0.0);
	return gasVoidFraction;
}
