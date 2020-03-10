//*****************************************************************************
// Author: Zack Taylor 
//
// Driver class for transported species. This class is used to build the 
// species transport problem
//*****************************************************************************
#ifndef SPECIESDRIVER_H
#define SPECIESDRIVER_H
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include "modelMesh.h"
#include "meshCellData.h"
#include "meshCellFace.h"
#include "species.h"
#include "matrixExponential.h"
#include "convectionLimiter.h"
#include "mpiProcess.h"
#include "matrixTypes.h"
#include "vectorTypes.h"
#include "exception.h"

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
	SparseMatrixD A;
	// Initial condition
	VectorD N0;
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
	// Sets the species concentration
	void setSpeciesCon(int, int, int, double);
	// Sets the species source terms
	void setSpeciesSource(int, int, int, std::vector<double>, double);
	// Sets a boundary condition in a cell
	void setBoundaryCondition(std::string, std::string, int, double = 0);
	// Call to make solver rebuild the A matrix before the next solve
	void resetMatrix();
	// Solves the transient species transport equation with matrix exp
	void solve(double);
	// Solves the transient species transport equation with implicit solve
	void solveImplicit(double);
	// Solves the steady state species transport equation
	void solve();
	// Gives ability to se the matrix exp solver
	void setMatrixExpSolver(std::string, bool = false, int = 10);
	// Cleans species
	void clean();

	private:
	// Exponential solver
	matrixExponential *expSolver;
	// Builds the transition matrix
	SparseMatrixD buildTransMatrix(bool, double);
	// Builds the initial condition vector
	VectorD buildInitialConditionVector(bool);
	// Builds the b vector (holding the constant sources) 
	VectorD buildbVector();
	// Unpacks the solution
	void unpackSolution(const VectorD&);
	// Gets the i or j index for transition matrix
	int getAi(int, int, int, int);
	// Calculates the species convection slope across in a cell
	double calcSpecConvectiveSlope(int, int, int, double);
	// Sets a generic boundary condition
	void setGeneralBoundaryCondition(std::string, int, int, double);
	// Sets a periodic boundary condiiton
	void setPeriodicBoundaryCondition(int);
};
#endif
