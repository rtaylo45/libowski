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
#include <istream>
#include <filesystem>
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
#include "ODEintegrator.h"
#include "massTransfer.h"

class speciesDriver {

	// Class attributes
	public:
	// Pointer to the model object
	modelMesh* modelPtr = nullptr;
	// Number of species in the model
	int numOfSpecs = 0;
	// Number of dummy species needs to be added. 
	int dummySpec = 1;
	// Logical set after the matrix has been built
	bool matrixInit = false;
	// Last solve time [s]
	double lastSolveTime = 0.0;
	// Transition matrix
	SparseMatrixD A;
	// Initial condition
	VectorD N0;
	// Convection flux limiter type
	fluxLimiter fluxLim = fluxLimiter(6);

	// Class methods
	public:
	// Constructor
	speciesDriver(modelMesh* model);
	// Adds a species to the system
	int addSpecies(double, double = 0.0, double = 0.0, std::string = "None",
		bool = true);
	// Adds species from a file generated from pyLibowski
	std::vector<int> addSpeciesFromFile(std::string);
	// Gets a pointer to the species object
	species* getSpeciesPtr(int, int, int);
	// Gets the species concentration
	double getSpecies(int, int, int);
	// Gets the species name
	std::string getSpeciesName(int, int , int);
	// Sets the species concentration
	void setSpeciesCon(int, int, int, double);
	// Sets the species source terms
	void setSpeciesSource(int, int, int, std::vector<double>, double = 0.0,
		std::vector<double> = std::vector<double>());
	// Sets the species source terms from files
	void setSpeciesSourceFromFile(std::string, std::string = "None");
	// Sets a boundary condition in a cell
	void setBoundaryCondition(std::string, std::string, int, double = 0);
	// Sets a boundary condition in a cell for a list of isotopes
	void setBoundaryCondition(std::string, std::string, std::vector<int>,
		std::vector<double> = std::vector<double>());
	// Call to make solver rebuild the A matrix before the next solve
	void resetMatrix();
	// Solves the transient species transport equation with matrix exp
	void solve(double);
	// Solves the transient species transport equation with implicit solve
	void solveImplicit(double);
	// Solves the steady state species transport equation
	void solve();
	// Gives ability to set the matrix exp solver
	void setMatrixExpSolver(std::string, bool = false, int = 10);
	// Gives ability to set the integrator solver
	void setIntegratorSolver(std::string, std::string);
	// Gives ability to set the flux limiter function
	void setFluxLimiter(std::string);
	// Writes out the base line transition matrix to a csv file. This matrix
	// is not multiplied by the time step size. This function must be called
	// After all of the speices source terms are set and all mesh parameters 
	// are set. Call this right before you call the solve method
	void writeTransitionMatrixToFile(std::string);
	// Wrties out the initial condition vector to a csv file
	void writeInitialConditionToFile(std::string);
	/// Sets the Krylov subspace dimension of the matexp solver
	void setKrylovSubspaceDimension(int dim);
	// Cleans species
	void clean();

	private:
	// Sets the decay coefficients
	void setDecaySource(int i, int j, int specID, std::string, std::vector<double>);
	// Sets the transmutation coefficients
	void setTransSource(int i, int j, int specID, std::string, std::vector<double>);
	// Solver step
	int step = 0;
	// Exponential solver
	matrixExponential *expSolver;
	// Integrator solver
	ODEintegrator *intSolver;
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
	double calcSpecConvectiveSlope(meshCell*, connection*, int, double);
	// Sets a generic boundary condition
	void setGeneralBoundaryCondition(std::string, int, int, double);
	// Sets a periodic boundary condiiton
	void setPeriodicBoundaryCondition(int);
	// Calculates the source term for deferred correction
	double calcDefCor(meshCell*, connection*, int, double);
	// Calculates a vector of source terms for the deferred corrections
	VectorD calcDefSourceVector();
};
#endif
