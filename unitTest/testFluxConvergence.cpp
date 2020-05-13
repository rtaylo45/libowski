#define _USE_MATH_DEFINES
#include <assert.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <cmath>

#include "speciesDriver.h"
#include "modelMesh.h"
#include "mpiProcess.h"
#include "matrixTypes.h"
#include "utilBase.h"


//*****************************************************************************
// All convergence test will use the following PDE but each test will have
// different initial conditions
//
// du/dt = -du/dx
// 
// Majority of test were taken from. 
//
// Ninth MSU-UAB Conference on Differential Equations and Computational Simulations.
// Electronic Journal of Differential Equations, Conference 20 (2013), pp. 53–63.
// ISSN: 1072-6691. URL: http://ejde.math.txstate.edu or http://ejde.math.unt.edu
// ftp ejde.math.txstate.edu
//
// TOTAL VARIATION STABILITY AND SECOND-ORDER ACCURACY AT EXTREMA
//*****************************************************************************


//*****************************************************************************
// Test one
//
// x		[0, 100]
// IC.	e^-((x-30)/10)^2
//*****************************************************************************
void problem1(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10, 20, 40, 80, 160, 320};
	//std::vector<int> numOfxCells{100, 200, 400, 800, 1600};
	double xLength = 100., yLength = 0.0;
	double dx;
	double tEnd = 2.0;
	double cSol, cCon;
	double t = 0.0;
	int numOfSteps = 200;
	double dt = tEnd/numOfSteps;
	double x1, x2, initCon, xc;
	double l1, l2, linf;	
	int cID;
	double velocity = 2.0;
	double maxRelativeError = 0.0, relativeError = 0.0, rmse = 0.0;
	double maxRelativeErrorOld = 0.0, rmseOld = 0.0, maxRelRate = 0.0, rmseRate = 0.0;
	std::string outputFileName, solverType = "hyperbolic";
	meshCell* cell = nullptr;
	// Sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "fluxProblem1.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "Total problem time: " << tEnd << "\n";
	outputFile << "Total problem length: " << xLength << "\n";
	outputFile << "Refinement: " << "space" << "\n";

	// loops over number of cells
	for (int &xCells : numOfxCells){
		//dx = xLength/(double)xCells;
		//double dt = 0.8*dx;
		//int numOfSteps = ceil(tEnd/dt);
		// Build the mesh
		modelMesh model = modelMesh(xCells, yCells, xLength, yLength);
		// Builds the species driver
		speciesDriver spec = speciesDriver(&model);
		// Add species
		cID = spec.addSpecies(1.0, 0.0, 0.0);
		// Sets the Matrix expential solver
		spec.setMatrixExpSolver(solverType);
		// x velocity
		model.setConstantXVelocity(velocity);
		// Boundary conditions
		spec.setBoundaryCondition("periodic", "west", cID);
		spec.setBoundaryCondition("periodic", "east", cID);
		outputFile << "Solver: " << solverType << "\n";
		outputFile << "dx: " << xLength/(double)xCells << "\n";
		outputFile << "dt: " << dt << "\n";
		outputFile << "variables " << "x " << "Solution " << "Libowski" << "\n";
		maxRelativeError = 0.0;
		rmse = 0.0;
		l1 = 0.0;
		l2 = 0.0;

		// Sets the intial condition and sources
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				cell = model.getCellByLoc(i,j);	

				// Calculates the x positions as the cell faces
				dx = cell->dx;
				xc = cell->x;
				x2 = xc + dx/2.;
				x1 = xc - dx/2.;

				// Calculates the initial concentration from MVT. 
				initCon = 1./dx*(5.*sqrt(M_PI)*(erf(3.-x1/10.) - erf(3.-x2/10.)));
				//initCon = exp(-std::pow((xc-30.)/10., 2.));

				spec.setSpeciesCon(i, j, cID, initCon);
			}
		}

		// Solves the problem for each time step
		for (int step = 1; step <= numOfSteps; step++){
			t = step*dt;
			spec.solve(t);
		}
		// Loops to get solution info
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				cell = model.getCellByLoc(i,j);	
				xc = cell->x;
				cSol = exp(-std::pow((xc-velocity*t-30.)/10., 2.));
				cCon = spec.getSpecies(i, j, cID);
				relativeError = std::abs(cSol-cCon);
				l1 += std::abs(cSol-cCon);
				l2 += l1*l1;
				maxRelativeError = std::max(maxRelativeError, relativeError);
				rmse += std::pow(relativeError,2.);

				outputFile << xc << " " << cSol << " " << cCon << "\n";
			}
		}
		printf("%15s %2.3f %4.2E %4.2e %4.2e %4.2e %4.2e \n", solverType.c_str(), dt, dx, 
				maxRelativeError, std::pow(rmse/float(xCells), 0.5), 
				l1/float(xCells), l2/float(xCells));
		outputFile << "\n";
	}
	outputFile << "end";
}

//*****************************************************************************
// Test 2
//
// x		[0, 100]
// IC.	0.0
//*****************************************************************************
void problem2(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10, 20, 40, 80, 160, 320};
	//std::vector<int> numOfxCells{100, 200, 400, 800, 1600};
	double xLength = 100., yLength = 0.0;
	double dx;
	double tEnd = 20.0;
	double cSol, cCon;
	double t = 0.0;
	int numOfSteps = 200;
	double dt = tEnd/numOfSteps;
	double x1, x2, initCon, xc;
	int cID;
	double velocity = 2.0;
	double maxRelativeError = 0.0, relativeError = 0.0, rmse = 0.0;
	std::string outputFileName, solverType = "hyperbolic";
	meshCell* cell = nullptr;
	// Sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "fluxProblem2.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "Total problem time: " << tEnd << "\n";
	outputFile << "Total problem length: " << xLength << "\n";
	outputFile << "Refinement: " << "space" << "\n";

	// loops over number of cells
	for (int &xCells : numOfxCells){
		//dx = xLength/(double)xCells;
		//double dt = 0.8*dx;
		//int numOfSteps = ceil(tEnd/dt);
		// Build the mesh
		modelMesh model = modelMesh(xCells, yCells, xLength, yLength);
		// Adds boundarys
		model.addBoundarySurface("east");
		model.addBoundarySurface("west");
		// Builds the species driver
		speciesDriver spec = speciesDriver(&model);
		// Add species
		cID = spec.addSpecies(1.0, 0.0, 0.0);
		// Sets the Matrix expential solver
		spec.setMatrixExpSolver(solverType);
		// x velocity
		model.setConstantXVelocity(velocity);
		// Boundary conditions
		spec.setBoundaryCondition("dirichlet", "west", cID, 1.0);
		spec.setBoundaryCondition("free flow", "east", cID);
		outputFile << "Solver: " << solverType << "\n";
		outputFile << "dx: " << xLength/(double)xCells << "\n";
		outputFile << "dt: " << dt << "\n";
		outputFile << "variables " << "x " << "Solution " << "Libowski" << "\n";
		maxRelativeError = 0.0;
		rmse = 0.0;

		// Solves the problem for each time step
		for (int step = 1; step <= numOfSteps; step++){
			t = step*dt;
			spec.solve(t);
		}
		// Loops to get solution info
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				cell = model.getCellByLoc(i,j);	
				xc = cell->x;
				if (xc < velocity*t){
					cSol = 1.0;
				}
				else {
					cSol = 0.0;
				}
				cCon = spec.getSpecies(i, j, cID);
				relativeError = std::abs(cSol-cCon);
				maxRelativeError = std::max(maxRelativeError, relativeError);
				rmse += std::pow(relativeError,2.);

				outputFile << xc << " " << cSol << " " << cCon << "\n";
			}
		}
		printf("%15s %2.3f %4.2E %4.2e %4.2e %4.2e %4.2e \n", solverType.c_str(), dt, dx, 
				maxRelativeError, std::pow(rmse/float(xCells), 0.5));
		outputFile << "\n";
	}
	outputFile << "end";
}

//*****************************************************************************
// Test 3
//
// x		[0, 100]
// IC.	0.0
//*****************************************************************************
void problem3(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{200};
	//std::vector<int> numOfxCells{100, 200, 400, 800, 1600};
	std::vector<std::string> fluxLimiters {"First order upwind", "superbee",
		"muscl"};
	double xLength = 100., yLength = 0.0;
	double dx;
	double tEnd = 5.0;
	double cSol, cCon;
	double t = 0.0;
	int numOfSteps = 200;
	double dt = tEnd/numOfSteps;
	double x1, x2, initCon, xc;
	int cID;
	double velocity = 2.0;
	double maxRelativeError = 0.0, relativeError = 0.0, rmse = 0.0;
	std::string outputFileName, solverType = "hyperbolic";
	meshCell* cell = nullptr;
	// Sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "fluxProblem3.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "Total problem time: " << tEnd << "\n";
	outputFile << "Total problem length: " << xLength << "\n";
	outputFile << "Refinement: " << "space" << "\n";

	for (std::string &limiterType : fluxLimiters){
		// loops over number of cells
		for (int &xCells : numOfxCells){
			// Build the mesh
			modelMesh model = modelMesh(xCells, yCells, xLength, yLength);
			// Adds boundarys
			model.addBoundarySurface("east");
			model.addBoundarySurface("west");
			// Builds the species driver
			speciesDriver spec = speciesDriver(&model);
			// Sets the limiter function
			spec.setFluxLimiter(limiterType);
			// Add species
			cID = spec.addSpecies(1.0, 0.0, 0.0);
			// Sets the Matrix expential solver
			spec.setMatrixExpSolver(solverType);
			// x velocity
			model.setConstantXVelocity(velocity);
			// Boundary conditions
			spec.setBoundaryCondition("periodic", "west", cID);
			spec.setBoundaryCondition("periodic", "east", cID);
			outputFile << "Solver: " << limiterType << "\n";
			outputFile << "dx: " << xLength/(double)xCells << "\n";
			outputFile << "dt: " << dt << "\n";
			outputFile << "variables " << "x " << "Solution " << "Libowski" << "\n";
			maxRelativeError = 0.0;
			rmse = 0.0;

			// Sets the intial condition and sources
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	

					// Calculates the x positions as the cell faces
					dx = cell->dx;
					xc = cell->x;
					x2 = xc + dx/2.;
					x1 = xc - dx/2.;

					if (xc < 20. or xc > 80.){
						initCon = 0.0;
					}
					else if (xc >= 20. and xc <= 40.){
						initCon = 1.0;
					}
					else if (xc >= 60. and xc <= 80.){
						initCon = 1.0;
					}
					else {
						initCon = 0.5;
					}


					spec.setSpeciesCon(i, j, cID, initCon);
				}
			}

			// Solves the problem for each time step
			for (int step = 1; step <= numOfSteps; step++){
				t = step*dt;
				spec.solve(t);
			}
			// Loops to get solution info
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	
					xc = cell->x;

					if (xc < 20.+velocity*t or xc > 80.+velocity*t){
						cSol = 0.0;
					}
					else if (xc >= 20.+velocity*t and xc <= 40.+velocity*t){
						cSol = 1.0;
					}
					else if (xc >= 60.+velocity*t and xc <= 80.+velocity*t){
						cSol = 1.0;
					}
					else {
						cSol = 0.5;
					}

					cCon = spec.getSpecies(i, j, cID);
					relativeError = std::abs(cSol-cCon);
					maxRelativeError = std::max(maxRelativeError, relativeError);
					rmse += std::pow(relativeError,2.);

					outputFile << xc << " " << cSol << " " << cCon << "\n";
				}
			}
			printf("%15s %2.3f %4.2E %4.2e %4.2e %4.2e %4.2e \n", solverType.c_str(), dt, dx, 
					maxRelativeError, std::pow(rmse/float(xCells), 0.5));
			outputFile << "\n";
		}
	}
	outputFile << "end";
}

//*****************************************************************************
// Test 4
//
// x		[0, 100]
// IC.	0.0
//*****************************************************************************
void problem4(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{200};
	//std::vector<int> numOfxCells{100, 200, 400, 800, 1600};
	std::vector<std::string> fluxLimiters {"First order upwind", "superbee"};
	double xLength = 100., yLength = 0.0;
	double dx;
	double tEnd = 10.0;
	double cSol, cCon;
	double t = 0.0;
	int numOfSteps = 200;
	double dt = tEnd/numOfSteps;
	double x1, x2, initCon, xc;
	int cID;
	double velocity = 2.0;
	double maxRelativeError = 0.0, relativeError = 0.0, rmse = 0.0;
	std::string outputFileName, solverType = "hyperbolic";
	meshCell* cell = nullptr;
	// Sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "fluxProblem4.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "Total problem time: " << tEnd << "\n";
	outputFile << "Total problem length: " << xLength << "\n";
	outputFile << "Refinement: " << "space" << "\n";

	for (std::string &limiterType : fluxLimiters){
		// loops over number of cells
		for (int &xCells : numOfxCells){
			// Build the mesh
			modelMesh model = modelMesh(xCells, yCells, xLength, yLength);
			// Adds boundarys
			model.addBoundarySurface("east");
			model.addBoundarySurface("west");
			// Builds the species driver
			speciesDriver spec = speciesDriver(&model);
			// Sets the limiter function
			spec.setFluxLimiter(limiterType);
			// Add species
			cID = spec.addSpecies(1.0, 0.0, 0.0);
			// Sets the Matrix expential solver
			spec.setMatrixExpSolver(solverType);
			// x velocity
			model.setConstantXVelocity(velocity);
			// Boundary conditions
			spec.setBoundaryCondition("periodic", "west", cID);
			spec.setBoundaryCondition("periodic", "east", cID);
			outputFile << "Solver: " << limiterType << "\n";
			outputFile << "dx: " << xLength/(double)xCells << "\n";
			outputFile << "dt: " << dt << "\n";
			outputFile << "variables " << "x " << "Solution " << "Libowski" << "\n";
			maxRelativeError = 0.0;
			rmse = 0.0;

			// Sets the intial condition and sources
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	

					// Calculates the x positions as the cell faces
					dx = cell->dx;
					xc = cell->x;
					x2 = xc + dx/2.;
					x1 = xc - dx/2.;

					if (xc < 20. or xc > 30.){
						initCon = 0.0;
					}
					else {
						initCon = 1.0;
					}

					spec.setSpeciesCon(i, j, cID, initCon);
				}
			}

			// Solves the problem for each time step
			for (int step = 1; step <= numOfSteps; step++){
				t = step*dt;
				spec.solve(t);
				//spec.solveImplicit(t);
			}
			// Loops to get solution info
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	
					xc = cell->x;

					if (xc < 20.+velocity*t or xc > 30.+velocity*t){
						cSol = 0.0;
					}
					else {
						cSol = 1.0;
					}

					cCon = spec.getSpecies(i, j, cID);
					relativeError = std::abs(cSol-cCon);
					maxRelativeError = std::max(maxRelativeError, relativeError);
					rmse += std::pow(relativeError,2.);

					outputFile << xc << " " << cSol << " " << cCon << "\n";
				}
			}
			printf("%15s %2.3f %4.2E %4.2e %4.2e %4.2e %4.2e \n", solverType.c_str(), dt, dx, 
					maxRelativeError, std::pow(rmse/float(xCells), 0.5));
			outputFile << "\n";
		}
	}
	outputFile << "end";
}
int main(){

	int myid = mpi.rank;
	int numprocs = mpi.size;

	problem1(myid);
	//problem2(myid);
	//problem3(myid);
	//problem4(myid);
}