//*****************************************************************************
// Author: Zack Taylor
//
// Surface type that holds species boundary conditions and allows for surface
// chemisty source or sink terms.
//*****************************************************************************
#ifndef SURFACE_H
#define SURFACE_H
#include "species.h"
#include <string>
#include <assert.h>

class surface {

	public:
	// Constructure
	surface();
	// Initilizes the surface
	void set();
	// Add species
	void addSpecies(double, double = 0.0, double = 0.0, std::string = "None", 
		bool = true);
	// Gets a pointer to the species
	species* getSpeciesPtr(int);
	// Clean
	void clean();
	// Logical to set if the surface is initilized (exist)
	bool isInit = false;

	private:
	// Vector of the species on the surface
	std::vector<species> speciesVector;

};
#endif
