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
//
// The units for these test are arbritrary as the problems are taken directly
// from papers
//*****************************************************************************


//*****************************************************************************
// Test one x direction
//
// v		2.0
// x		[0, 100]
// I.C.	e^-((x-30)/10)^2
// B.C.	periodic
//*****************************************************************************
void problem1x(int myid, std::string solverType){
	int yCells = 1;
	std::vector<int> numOfxCells{10, 20, 40, 80, 160, 320};
	std::vector<double> linfVector, l1Vector, l2Vector;
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
	double absError = 0.0;
	double linfRate, l2Rate, l1Rate;
	meshCell* cell = nullptr;
	// Sets the ouput file name
	std::ofstream outputFile;
	std::string outputFileName = "fluxProblem1x.out";
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
		linf = 0.0;
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
		if (myid == 0){
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	
					xc = cell->x;
					cSol = exp(-std::pow((xc-velocity*t-30.)/10., 2.));
					cCon = spec.getSpecies(i, j, cID);
					absError = std::abs(cSol-cCon);
					l1 += std::abs(cSol-cCon);
					l2 += l1*l1;
					linf = std::max(linf, absError);

					outputFile << xc << " " << cSol << " " << cCon << "\n";
				}
			}
			//printf("%15s %2.3f %4.2E %4.2e %4.2e %4.2e \n", solverType.c_str(), dt, dx, 
					//linf, l1/float(xCells), l2/float(xCells));
			outputFile << "\n";
			linfVector.push_back(linf);
			l1Vector.push_back(l1/float(xCells));
			l2Vector.push_back(l2/float(xCells));
		}
	}
	outputFile << "end";
	// Loop through error vectors to check convergence
	if (myid == 0){
		for (int i = 1; i < linfVector.size(); i++){
			linfRate = log2(linfVector[i-1]/linfVector[i]);
			l1Rate = log2(l1Vector[i-1]/l1Vector[i]);
			l2Rate = log2(l2Vector[i-1]/l2Vector[i]);
			int xCells = numOfxCells[i];
			printf("%15s %4d %4.2f %4.2f %4.2f \n", solverType.c_str(), xCells, 
				linfRate, l1Rate, l2Rate);
			l1Rate = round(l1Rate);
			l2Rate = round(l2Rate);
			assert(l1Rate == 2);
			assert(l2Rate == 2);
		}
	}
}

//*****************************************************************************
// Test one y direction
//
// v		2.0
// x		[0, 100]
// IC.	e^-((x-30)/10)^2
// B.C.	periodic
//*****************************************************************************
void problem1y(int myid, std::string solverType){
	int xCells = 1;
	std::vector<int> numOfyCells{10, 20, 40, 80, 160, 320};
	std::vector<double> linfVector, l1Vector, l2Vector;
	double xLength = 0., yLength = 100.0;
	double dy;
	double tEnd = 2.0;
	double cSol, cCon;
	double t = 0.0;
	int numOfSteps = 200;
	double dt = tEnd/numOfSteps;
	double y1, y2, initCon, yc;
	double l1, l2, linf;	
	int cID;
	double velocity = 2.0;
	double absError = 0.0;
	double linfRate, l2Rate, l1Rate;
	meshCell* cell = nullptr;
	// Sets the ouput file name
	std::ofstream outputFile;
	std::string outputFileName = "fluxProblem1y.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "Total problem time: " << tEnd << "\n";
	outputFile << "Total problem length: " << yLength << "\n";
	outputFile << "Refinement: " << "space" << "\n";
	
	// loops over number of cells
	for (int &yCells : numOfyCells){
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
		model.setConstantYVelocity(velocity);
		// Boundary conditions
		spec.setBoundaryCondition("periodic", "north", cID);
		spec.setBoundaryCondition("periodic", "south", cID);
		outputFile << "Solver: " << solverType << "\n";
		outputFile << "dx: " << yLength/(double)yCells << "\n";
		outputFile << "dt: " << dt << "\n";
		outputFile << "variables " << "x " << "Solution " << "Libowski" << "\n";
		linf = 0.0;
		l1 = 0.0;
		l2 = 0.0;

		// Sets the intial condition and sources
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				cell = model.getCellByLoc(i,j);	

				// Calculates the x positions as the cell faces
				dy = cell->dy;
				yc = cell->y;
				y2 = yc + dy/2.;
				y1 = yc - dy/2.;

				// Calculates the initial concentration from MVT. 
				initCon = 1./dy*(5.*sqrt(M_PI)*(erf(3.-y1/10.) - erf(3.-y2/10.)));

				spec.setSpeciesCon(i, j, cID, initCon);
			}
		}

		// Solves the problem for each time step
		for (int step = 1; step <= numOfSteps; step++){
			t = step*dt;
			spec.solve(t);
		}
		// Loops to get solution info
		if (myid == 0){
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	
					yc = cell->y;
					cSol = exp(-std::pow((yc-velocity*t-30.)/10., 2.));
					cCon = spec.getSpecies(i, j, cID);
					absError = std::abs(cSol-cCon);
					l1 += std::abs(cSol-cCon);
					l2 += l1*l1;
					linf = std::max(linf, absError);

					outputFile << yc << " " << cSol << " " << cCon << "\n";
				}
			}
			//printf("%15s %2.3f %4.2E %4.2e %4.2e %4.2e \n", solverType.c_str(), dt, dy, 
					//linf, l1/float(yCells), l2/float(yCells));
			outputFile << "\n";
			linfVector.push_back(linf);
			l1Vector.push_back(l1/float(yCells));
			l2Vector.push_back(l2/float(yCells));
		}
	}
	outputFile << "end";
	// Loop through error vectors to check convergence
	if (myid == 0){
		for (int i = 1; i < linfVector.size(); i++){
			linfRate = log2(linfVector[i-1]/linfVector[i]);
			l1Rate = log2(l1Vector[i-1]/l1Vector[i]);
			l2Rate = log2(l2Vector[i-1]/l2Vector[i]);
			int yCells = numOfyCells[i];
			printf("%15s %4d %4.2f %4.2f %4.2f \n", solverType.c_str(), yCells, 
				linfRate, l1Rate, l2Rate);
			l1Rate = round(l1Rate);
			l2Rate = round(l2Rate);
			assert(l1Rate == l1Rate);
			assert(l2Rate == l2Rate);
		}
	}
}
//*****************************************************************************
// Test 2
//
// v		2.0
// x		[0, 100]
// I.C.	0.0
// B.C.  west: Dirichlet 1.0
//			east: Free flow
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
		if (myid == 0){
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
			printf("%15s %2.3f %4.2E %4.2e %4.2e \n", solverType.c_str(), dt,
			dx, maxRelativeError, std::pow(rmse/float(xCells), 0.5));
			outputFile << "\n";
		}
	}
	outputFile << "end";
}

//*****************************************************************************
// Test 3
//
// v		2.0
// x		[0, 100]
// IC.	1.0 for 20 <= x <= 40
//			1.0 for 60 <= x <= 80
//			0.5 for 40 <  x <  60
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
			if (myid == 0){
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
				printf("%15s %2.3f %4.2E %4.2e %4.2e \n", solverType.c_str(),
				dt, dx, maxRelativeError, std::pow(rmse/float(xCells), 0.5));
				outputFile << "\n";
			}
		}
	}
	outputFile << "end";
}

//*****************************************************************************
// Test 4
//
// v		2.0
// x		[0, 100]
// I.C.	1.0 for 20 <= x <= 30
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
			if (myid == 0){
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
				printf("%15s %2.3f %4.2E %4.2e %4.2e \n", solverType.c_str(),
				dt, dx, maxRelativeError, std::pow(rmse/float(xCells), 0.5));
				outputFile << "\n";
			}
		}
	}
	outputFile << "end";
}

//*****************************************************************************
// 2 species 1D diffusion problem, taken from the following paper
// NUMERICAL METHODS FOR STIFF REACTION-DIFFUSION SYSTEMS
//	By: Chou, Zhang, Zhao, and Nie
//
//	Diff eqs:
//		dU/dt = d*Uxx - a*U + V
//		dV/dt = d*Vxx - b*V
//
//	Solution:
//		U(x,t) = (e^(-(a+d)*t) + e^(-(b+d)*t))*cos(x).
//		U(x,t) = (a-b)*e^(-(b+d)*t)*cos(x).
//
// They don't explicitly give the initial condiiton in the paper, but if you
// plug in zero to the solution you can get it for U and V.
//
//	BC's:
//		Ux(0,t) = 0, Uv(0.t) = 0, U(pi/2,t) = 0, V(pi/2,t) = 0
//
//	Note: The units for this problem kinda dont matter either. This is just an
//	internal ODE test, the ODE solution is provided and calculated using the 
//	same problem units.
//*****************************************************************************
void problem5(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10, 20, 40, 80, 160, 320};
	double xLength = M_PI/2., yLength = 0.0;
	double numOfSteps = 1;
	double tEnd = 1.0;
	double dt = tEnd/numOfSteps, t = 0;
	double a = 0.1, b = 0.01, d = 1.0;	// Diffusion dominated
	//double a = 2.0, b = 1.0, d = 0.001;		// Reaction dominated
	//double a = 100.0, b = 1.0, d = 0.001;		// Stiff reaction dominated
	double UCon, VCon, USol, VSol;
	int UID, VID;
	double x1, x2, xc, initCon, x, dx;
	double l1, l2, linf;
	double linfRate, l2Rate, l1Rate;
	double absError, runtime;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> Ucoeffs = {-a, 1.0};
	std::vector<double> Vcoeffs = {0.0, -b};
	std::vector<std::string> solvers {"CRAM", "parabolic", "hyperbolic", 
		"pade-method1", "pade-method2", "taylor"};

	// Loops over different solvers
	for (std::string &solverType : solvers){
		std::vector<double> linfVector, l1Vector, l2Vector, runTimeVector;

		std::ofstream outputFile;
		outputFileName = "problem5"+solverType+".out";
		outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
		
		// Loops over number of cells
		for (int &xCells : numOfxCells){

			// Build the Mesh
			modelMesh model(xCells, yCells, xLength, yLength);
			// Sets the surfaces
			model.addBoundarySurface("east");
			model.addBoundarySurface("west");
			// Build species driver
			speciesDriver spec = speciesDriver(&model);

			// Add species. I will add the initial condition later
			UID = spec.addSpecies(1.0, 0.0, d);
			VID = spec.addSpecies(1.0, 0.0, d);
			// Sets the species matrix exp solver
			spec.setMatrixExpSolver(solverType);

			// Add BCs
			spec.setBoundaryCondition("newmann","west", UID, 0.0);
			spec.setBoundaryCondition("dirichlet","east", UID, 0.0);
			spec.setBoundaryCondition("newmann","west", VID, 0.0);
			spec.setBoundaryCondition("dirichlet","east", VID, 0.0);

			// Sets the intial condition
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	

					// Calculates the x positions as the cell faces
					dx = cell->dx;
					xc = cell->x;
					x2 = xc + dx/2;
					x1 = xc - dx/2;

					// Calculates the initial concentration from MVT. 
					initCon = (1./dx)*(sin(x2) - sin(x1));		

					spec.setSpeciesCon(i,j,UID, 2*initCon);
					spec.setSpeciesCon(i,j,VID, (a-b)*initCon);

					// Sets the sourses
					spec.setSpeciesSource(i, j, UID, Ucoeffs);
					spec.setSpeciesSource(i, j, VID, Vcoeffs);
				}
			}

			for (int step = 1; step <= numOfSteps; step++){
				t = step*dt;
				// Solve with CRAM
				auto start = std::chrono::high_resolution_clock::now();
				spec.solve(t);
				auto end = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
					end - start);

				linf = 0.0; l1 = 0.0; l2 = 0.0;
				// Gets species Concentrations
				if (myid==0){
					for (int i = 0; i < xCells; i++){
						for (int j = 0; j < yCells; j++){
							cell = model.getCellByLoc(i,j);	

							// Caclulate analytical solution
							x = cell->x;
							USol = (exp(-(a+d)*t) + exp(-(b+d)*t))*cos(x);
							VSol = (a-b)*exp(-(b+d)*t)*cos(x);
							// Get libowski solution
							UCon = spec.getSpecies(i, j, UID);
							VCon = spec.getSpecies(i, j, VID);
						
							absError = std::abs(USol-UCon)/USol + std::abs(VSol-VCon)/VSol;
							l1 += absError;
							l2 += l1*l1;
							linf = std::max(linf, absError);

						}
					}
				}
				linfVector.push_back(linf);
				l1Vector.push_back(l1/float(xCells));
				l2Vector.push_back(l2/float(xCells));
				runTimeVector.push_back(duration.count()/1.e6);
			}
			model.clean();
			spec.clean();
		}
		// Loop through error vectors to check convergence
		if (myid == 0){
			for (int i = 1; i < linfVector.size(); i++){
				runtime = runTimeVector[i];	
				linfRate = log2(linfVector[i-1]/linfVector[i]);
				l1Rate = log2(l1Vector[i-1]/l1Vector[i]);
				l2Rate = log2(l2Vector[i-1]/l2Vector[i]);
				int xCells = numOfxCells[i];
				printf("%15s %4d %4.2f %4.2f %4.2f %4.2e %4.2e %4.2e %4.2e \n",
					solverType.c_str(), xCells, linfRate, l1Rate, l2Rate, 
					linfVector[i], l1Vector[i], l2Vector[i], runtime);
				l1Rate = round(l1Rate);
				l2Rate = round(l2Rate);
				assert(l1Rate == 2);
				assert(l2Rate == 2);
			}
		}
	}

}
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;
	std::vector<std::string> solvers {"CRAM", "parabolic", "hyperbolic", 
		"pade-method1", "pade-method2", "taylor"};

	// Loops over different solvers
	for (std::string &solverType : solvers){
		problem1x(myid, solverType);
		problem1y(myid, solverType);
	}
	problem2(myid); problem3(myid); problem4(myid); problem5(myid);

	mpi.finalize();
}
