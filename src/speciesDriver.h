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
#include <iostream>
#include <algorithm>
#include "modelMesh.h"
#include "meshCellData.h"
#include "meshCellFace.h"
#include "species.h"
#include "CRAM.h"
#include "convectionLimiter.h"
#include "mpiProcess.h"

class speciesDriver {

	// Class attributes
	public:
	// Pointer to the model object
	modelMesh* modelPtr = nullptr;
	// Number of species in the model
	int numOfSpecs = 0;
	// Number of dummy species needs to be added. 
	int dummySpec = 0;
	// Logical set after the matrix has been built
	bool matrixInit = false;
	// Last solve time 
	double lastSolveTime = 0.0;
	// Transition matrix
	Eigen::SparseMatrix<double> A;
	// Initial condition
	Eigen::VectorXd N0;
	// Convection flux limiter type
	fluxLimiter fluxLim = fluxLimiter(0);

	// Class methods
	public:
	// Constructor
	speciesDriver(modelMesh* model);
	// Adds a species to the system
	int addSpecies(double, double, double);
	// Gets a pointer to the species object
	species* getSpeciesPtr(int, int, int);
	// Gets the species concentration
	double getSpecies(int, int, int);
	// Sets the species source terms
	void setSpeciesSource(int, int, int, std::vector<double>, double);
	// Sets a boundary condition in a cell
	void setBoundaryCondition(int, int, int, double);
	// Call to make solver rebuild the A matrix before the next solve
	void resetMatrix();
	// Solves the species transport equation
	void solve(double);
	// Cleans species
	void clean();

	private:
	// Builds the transition matrix
	Eigen::SparseMatrix<double> buildTransMatrix();
	// Builds the initial condition vector
	Eigen::VectorXd buildInitialConditionVector();
	// Unpacks the solution
	void unpackSolution(Eigen::VectorXd);
	// Gets the i or j index for transition matrix
	int getAi(int, int, int, int);
	// Calculates the species convection slope across in a cell
	double calcSpecConvectiveSlope(int, int, int, double);
};
#endif
