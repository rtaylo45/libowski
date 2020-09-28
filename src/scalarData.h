//*****************************************************************************
// Author: Zack Taylor
//
//	Data object that holds scalar data in a mesh cell (except for species).
//*****************************************************************************
#ifndef SCALAR_H
#define SCALAR_H
#include <cassert>

class scalarData {

	// Member attributes
	private:
	// Temperature of cell in kelvin
	double T = -1.;
	// Presure in Pa
	double P = -1.;
	// Scalar neutron flux 1/cm^2/s
	double phi  = -1.;
	// Interfacial area concentration [1/m]
	double intAreaCon = -1.;

	// Member functions
	public:
	// Constructor
	scalarData();
	// Sets the temperature
	void setTemperature(double);
	// Sets the pressure
	void setPressure(double);
	// Sets the scalar neutron flux
	void setNeutronFlux(double);
	// Sets the interfacial area concentration 
	void setInterfacialAreaCon(double);
	// Gets the temperature
	double getTemperature();
	// Gets the pressure
	double getPressure();
	// Gets the scalar neutron flux
	double getNeutronFlux();
	// Gets the interfacial area concentration 
	double getInterfacialAreaCon();

};

#endif
