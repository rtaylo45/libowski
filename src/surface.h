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
	// Add species
	void addSpecies(double, double, double);
	// Gets a pointer to the species
	species* getSpeciesPtr(int);
	std::string name = "zack";
	// Clean
	void clean();

	private:
	// Vector of the species on the surface
	std::vector<species> speciesVector;

};
#endif
