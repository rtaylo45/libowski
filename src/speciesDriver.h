//*****************************************************************************
// Author: Zack Taylor 
//
// Driver class for transported species. This class is used to build the 
// species transport problem
//*****************************************************************************
#ifndef SPECIESDRIVER_H
#define SPECIESDRIVER_H
#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "modelMesh.h"
#include "meshCellData.h"
#include "meshCellFace.h"
#include "species.h"
#include "CRAM.h"

class speciesDriver {

	// Class attributes
	public:
	// Pointer to the model object
	modelMesh* modelPtr = nullptr;
	// Number of species in the model
	int numOfSpecs = 0;
	// Number of dummy species needs to be added. 
	int dummySpec = 0;

	// Class methods
	public:
	// Constructor
	speciesDriver(modelMesh* model);
	// Adds a species to the system
	int addSpecies(double, double);
	// Gets a pointer to the species object
	species* getSpeciesPtr(int, int, int);
	// Gets the species concentration
	double getSpecies(int, int, int);
	// Sets the species source terms
	void setSpeciesSource(int, int, int, std::vector<double>, double);
	// Solves the species transport equation
	void solve(double, Eigen::SparseVector<double>);
	// Cleans species
	void clean();

	private:
	// Builds the transition matrix
	Eigen::SparseMatrix<double> buildTransMatrix();
};
#endif
