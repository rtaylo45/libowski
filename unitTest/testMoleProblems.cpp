#define _USE_MATH_DEFINES
#include <string>
#include <math.h>
#include <cmath>
#include <iostream>
#include <chrono>
#include <stdio.h>

#include "mpiProcess.h"
#include "modelMesh.h"
#include "speciesDriver.h"
#include "meshCellData.h"
#include "species.h"
#include "utilBase.h"

//*****************************************************************************
// Mole problem 1
//
// This problem simple models a uniform distribution of a contrived species 
// that decays with a constant rate. The physical layout is a 1D model of 
// 100 cm pipe. the species does not move, resulting in a solution that is 
// only time-dependent.
//
//	Problem equations:
//			dC/dt = -lambda*C(x,t)
//
//	Domaine:
//			x = [0, 100]	cm
//			t = [0, 20]		s
//			lambda = 0.01	1/s
//
//	Initial conditions and BC's:
//			C(x, 0) = 1000.(x + 1)	kg/m^3
//			C(0,t) = C(100,t)			kg/m^3
//
//	Solution:
//			C(x,t) = C(x,0)*e^(-lambda*t)
//
//*****************************************************************************
void moleProblem1(int myid){
	int yCells = 1, xCells = 1000;
	std::vector<double> steps = {1, 2, 4, 8, 20, 40, 80, 200, 400};
	std::vector<std::string> solvers {"hyperbolic", "pade-method2", "taylor"};
	// convert cm to m
	double xLength = 100./100.; // m
	double yLength = 0.0; // m
	double tEnd = 20.0;	// seconds
	double lambda = 0.1;	// 1/s
	double cSol, cCon;
	int cID;
	double x1, x2, initCon, dx, xc, dt, t;
	double maxRelativeError = 0.0, relativeError = 0.0, rmse = 0.0;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> cCoeffs = {-lambda};
	// Sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "moleProblem1.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "Total problem time: " << tEnd << "\n";

	// Build the Mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	// Build species driver
	speciesDriver spec = speciesDriver(&model);

	// Loops over different solvers
	for (std::string &solverType : solvers){
		// Sets the species matrix exp solver
		spec.setMatrixExpSolver(solverType);


		// Loops over number of time steps
		for (double &numOfSteps	: steps){
			maxRelativeError = 0.0;
			rmse = 0.0;
			dt = tEnd/numOfSteps;
			outputFile << "Solver: " << solverType << "\n";
			outputFile << "dt: " << dt << "\n";
			outputFile << "x"	<< " " << "C" << "\n";

			// Add specs
			cID = spec.addSpecies(1.0, 0.0, 0.0);

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
					initCon = 1000.*((x2 + x1)/2. + 1.);
					double s = 1000.*(xc+1);

					spec.setSpeciesCon(i, j, cID, initCon);

					// Sets the sources
					spec.setSpeciesSource(i, j, cID, cCoeffs);
				}
			}
			t = 0.0;
			// Solve the problem
			auto start = std::chrono::high_resolution_clock::now();
			for (int step = 1; step <= numOfSteps; step++){
				t = step*dt;
				spec.solve(t);
			}
			auto end = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
					end - start);

			// Gets species Concentrations
			if (myid==0){
				for (int i = 0; i < xCells; i++){
					for (int j = 0; j < yCells; j++){
						cell = model.getCellByLoc(i,j);	

						// Caclulate analytical solution
						xc = cell->x;
						cSol = 1000.*(xc+1.)*exp(-lambda*t);

						// Get libowski solution
						cCon = spec.getSpecies(i, j, cID);

						relativeError = std::abs(cSol-cCon)/cSol;
						maxRelativeError = std::max(maxRelativeError, relativeError);
						rmse += std::pow(relativeError,2.);
						// xc is converted from m to cm
						outputFile << xc*100. << " " << cCon << "\n";
						assert(isApprox(cCon, cSol));
					}
				}
				outputFile << "\n";
				printf("%15s %4.2f %4.2e %4.2e %3.5f \n", solverType.c_str(), dt, 
					relativeError, std::pow(rmse/float(xCells), 0.5), duration.count()/1.e6);
			}
			// Clean species
			spec.clean();
		}
		spec.clean();
	}
}

//*****************************************************************************
// Mole problem 2
//
// This problem is the same as Problem 1, but with a constant flow rate of 
// 2 cm/s in the positive x direction.
//
//	Problem equations:
//			dC/dt = -vdC/dx -lambda*C(x,t)
//
//	Domaine:
//			x = [0, 100]	cm
//			t = [0, 20]		s
//			v = 2.0			cm/s
//			lambda = 0.01	1/s
//
//	Initial conditions and BC's:
//			C(x, 0) = 1000.0	kg/m^3
//			C(0,t) = 1000.0	kg/m^3
//
//	Solution:
//			C(x,t) = TBD
//
//*****************************************************************************
void moleProblem2(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10, 100, 1000};
	std::vector<double> steps = {1, 2, 4, 8, 20, 40, 80, 200, 400};
	std::vector<std::string> solvers {"hyperbolic","pade-method2", "taylor"};
	double xLength = 100./100.; // m
	double yLength = 0.0; // m
	double tEnd = 20.0;	// seconds
	double lambda = 0.01;	// 1/s
	double velocity = 2.0/100.; // m/s
	double cSol, cCon;
	int cID;
	double x1, x2, initCon, dx, xc, dt, t;
	double maxRelativeError = 0.0, relativeError = 0.0, rmse = 0.0;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> ccoeffs = {-lambda};
	// sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "moleproblem2.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "Total problem time: " << tEnd << "\n";
	outputFile << "Total problem length: " << xLength << "\n";
	outputFile << "Refinement: " << "time" << "\n";

	// Loops over different solvers
	for (std::string &solverType : solvers){

		// loops over number of cells
		for (int &xCells : numOfxCells){
			// build the mesh
			modelMesh model(xCells, yCells, xLength, yLength);
			// Add BC surface
			model.addBoundarySurface("west");
			model.addBoundarySurface("east");
			// build species driver
			speciesDriver spec = speciesDriver(&model);
			// set x velocity
			model.setConstantXVelocity(velocity);
			// Sets the species matrix exp solver
			spec.setMatrixExpSolver(solverType);

			// loops over number of time steps
			for (double &numofsteps	: steps){
				maxRelativeError = 0.0;
				rmse = 0.0;
				dt = tEnd/numofsteps;
				outputFile << "Solver: " << solverType << "\n";
				outputFile << "dx: " << xLength/(double)xCells << "\n";
				outputFile << "dt: " << dt << "\n";
				outputFile << "variables " << "x " << "Solution " << "Libowski" << "\n";

				// add specs
				cID = spec.addSpecies(1.0, 0.0, 0.0);

				spec.setBoundaryCondition("dirichlet","west", cID, 1000.0);
				spec.setBoundaryCondition("newmann","east", cID, 0.0);
				//spec.setBoundaryCondition("free flow","east", cID);

				// sets the intial condition and sources
				for (int i = 0; i < xCells; i++){
					for (int j = 0; j < yCells; j++){
						// calculates the initial concentration
						spec.setSpeciesCon(i, j, cID, 1000.0);
						// sets the sources
						spec.setSpeciesSource(i, j, cID, ccoeffs);
					}
				}
				t = 0.0;
				// solve the problem
				auto start = std::chrono::high_resolution_clock::now();
				for (int step = 1; step <= numofsteps; step++){
					t = step*dt;
					// solve with cram
					spec.solve(t);
				}
				auto end = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
						end - start);
				// gets species concentrations
				if (myid==0){
					for (int i = 0; i < xCells; i++){
						for (int j = 0; j < yCells; j++){
							cell = model.getCellByLoc(i,j);	

							// caclulate analytical solution
							xc = cell->x;
							if (xc < velocity*t){
								cSol = 1000.*exp(-lambda*(xc)/velocity);
							}
							else{
								cSol = 1000.*exp(-lambda*t);
							}

							// get libowski solution
							cCon = spec.getSpecies(i, j, cID);

							if (xCells == 1000 and numofsteps == 400.){
								assert(isApprox(cSol, cCon, 1.0, 1e-3));
							}
							//linfError = std::max(linfError, std::abs(cSol-cCon)/cSol);
							outputFile << xc << " " << cSol << " " << cCon << "\n";
							//std::cout << xc << std::endl; 
							// << cSol << " " << cCon << std::endl;
							relativeError = std::abs(cSol-cCon)/cSol;
							maxRelativeError = std::max(maxRelativeError, relativeError);
							rmse += std::pow(relativeError,2.);
						}
					}
					outputFile << "\n";
					//std::cout << solverType << " " << dx << " " << dt 
					//	<< " " << percentError << "\n";
					dx = xLength/(double)xCells;
					printf("%15s %2.3f %4.2E %4.2e %4.2e %3.5f \n", solverType.c_str(), 
					dt, dx, maxRelativeError, std::pow(rmse/float(xCells), 
					0.5), duration.count()/1.e6);
				}
				spec.clean();
			}
			spec.clean();
			model.clean();
		}
	}
	outputFile << "end";
}

//*****************************************************************************
// Mole problem 3
//
// This problem is the same as Problem 1, expect that instead of radioactive 
// decay, an Arrhenius deposition rate is applied instead. All removal from
// the liquid deposits on the wall and remains. None of the material is 
// deposited on the wall at the beginning.
//
//	Problem equations:
//			dCi/dt = -v*dCi/dx - lambda*Ci
//			dCw/dt = lambda*Cw
//
//	Domaine:
//			x = [0, 100]	cm 
//			t = [0, 20]		s
//			lambda = 0.1	1/s
//			v = 2.0			cm/s
//
//	Initial conditions and BC's:
//			Ci(x, 0) = 1000.0
//			Ci(0,t) = 1000.0
//			Cw(x, 0) = 0
//			Cw(0,t) = 0
//
//	Solution:
//			Ci(x,t) = 1000*e^(-lambda*t),			x >= v*t
//					  = 1000*e^(-lambda*x/v),		x < v*t
//
//			Cw(x,t) = 1000*(1 - e^(-lambda*t)),	x >= v*t
//					  = 1000*(1 - e^(-lambda*x/v) + lambda*e^(-lambda*x/v)*(t - x/v))
//						 x < v*t
//
//*****************************************************************************
void moleProblem3(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10, 100, 1000};
	//std::vector<double> steps = {400};
	std::vector<double> steps = {1, 2, 4, 8, 20, 40, 80, 200, 400};
	std::vector<std::string> solvers {"hyperbolic","pade-method2", "taylor"};
	double xLength = 100./100.; // m
	double yLength = 0.0; // m
	double tEnd = 20.0;	// seconds
	double lambda = 0.1;	// 1/s
	double velocity = 2.0/100.; // m/s
	double clSol, cwSol, clCon, cwCon;
	int clID, cwID;
	double x1, x2, initCon, dx, xc, dt, t;
	double maxRelativeError = 0.0, relativeError = 0.0, rmse = 0.0;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> clcoeffs = {-lambda, 0.0};
	std::vector<double> cwcoeffs = {lambda, 0.0};
	// sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "moleproblem3.out";
	if (myid == 0){
		outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
		outputFile << "Total problem time: " << tEnd << "\n";
		outputFile << "Total problem length: " << xLength << "\n";
		outputFile << "Refinement: " << "time" << "\n";
	}

	// Loops over different solvers
	for (std::string &solverType : solvers){
		// loops over number of cells
		for (int &xCells : numOfxCells){
			// build the mesh
			modelMesh model(xCells, yCells, xLength, yLength);
			// Add BC surface
			model.addBoundarySurface("west");
			model.addBoundarySurface("east");
			// build species driver
			speciesDriver spec = speciesDriver(&model);
			// set x velocity
			model.setConstantXVelocity(velocity);
			// Sets the species matrix exp solver
			spec.setMatrixExpSolver(solverType);

			// loops over number of time steps
			for (double &numofsteps	: steps){
				maxRelativeError = 0.0;
				rmse = 0.0;
				dt = tEnd/numofsteps;
				if (myid == 0){
					outputFile << "Solver: " << solverType << "\n";
					outputFile << "dx: " << xLength/(double)xCells << "\n";
					outputFile << "dt: " << dt << "\n";
					outputFile << "variables " << "x " << "clSol " << "clLib " << 
						"cwSol " << "cwLib" << "\n";
				}
				// add specs
				clID = spec.addSpecies(1.0, 0.0, 0.0, "cliquid", true);
				cwID = spec.addSpecies(1.0, 0.0, 0.0, "cwall", false);
				// Set BCs
				spec.setBoundaryCondition("dirichlet","west", clID, 1000.0);
				spec.setBoundaryCondition("dirichlet","west", cwID, 0.0);
				spec.setBoundaryCondition("newmann","east", clID, 0.0);
				spec.setBoundaryCondition("newmann","east", cwID, 0.0);

				//spec.setBoundaryCondition("free flow","east", clID);
				//spec.setBoundaryCondition("free flow","east", cwID);

				// sets the intial condition and sources
				for (int i = 0; i < xCells; i++){
					for (int j = 0; j < yCells; j++){
						// calculates the initial concentration
						spec.setSpeciesCon(i, j, clID, 1000.0);
						// sets the sources
						spec.setSpeciesSource(i, j, clID, clcoeffs);
						spec.setSpeciesSource(i, j, cwID, cwcoeffs);
					}
				}
				t = 0.0;
				auto start = std::chrono::high_resolution_clock::now();
				for (int step = 1; step <= numofsteps; step++){
					t = step*dt;
					spec.solve(t);
				}
				auto end = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
						end - start);
				// gets species concentrations
				if (myid==0){
					for (int i = 0; i < xCells; i++){
						for (int j = 0; j < yCells; j++){
							cell = model.getCellByLoc(i,j);	

							// caclulate analytical solution
							xc = cell->x;
							dx = cell->dx;
							if (xc < velocity*t){
								clSol = 1000.0*exp(-lambda*xc/velocity);
								cwSol = 1000.0*(1. + (lambda*(t-xc/velocity) - 
									1.)*exp(-lambda*xc/velocity));
							}
							else{
								clSol = 1000.0*exp(-lambda*t);
								cwSol = 1000.0*(1.-exp(-lambda*t));
							}
							// get libowski solution
							clCon = spec.getSpecies(i, j, clID);
							cwCon = spec.getSpecies(i, j, cwID);

							if (xCells == 1000 and numofsteps == 400.){
								assert(isApprox(cwSol, cwCon, 5.0, 1e-2));
								assert(isApprox(clSol, clCon, 5.0, 1e-2));
							}
							outputFile << xc << " " << clSol << " " << clCon <<
								" " << cwSol << " " << cwCon << "\n";
							relativeError = std::abs(clSol-clCon)/clSol + 
								std::abs(cwSol-cwCon)/cwSol;
							maxRelativeError = std::max(maxRelativeError, relativeError);
							rmse += std::pow(relativeError,2.);
						}
					}
					outputFile << "\n";
					printf("%15s %2.3f %4.2E %4.2e %4.2e %3.5f \n", 
						solverType.c_str(), dt, dx,maxRelativeError, 
						std::pow(rmse/float(2*xCells), 0.5), duration.count()/1.e6);
				}
				spec.clean();
			}
			spec.clean();
			model.clean();
		}
	}
	outputFile << "end";

}

//*****************************************************************************
// Mole problem 4
//
// This problem models the drift and decay of the 6 delayed neutron precursor 
// groups.
//
//	Problem equations:
//
//	Domaine:
//
//	Initial conditions and BC's:
//
//	Solution: 
//
//*****************************************************************************
void moleProblem4(int myid){
}

//*****************************************************************************
// Mole problem 5
//
// This problem is the same as Problem 3, but with the liquid flowing at a 
// constant velocity in the positive x direction.
//
//	Problem equations:
//
//	Domaine:
//
//	Initial conditions and BC's:
//
//	Solution: 
//
//*****************************************************************************
void moleProblem5(int myid){
}

//*****************************************************************************
// Mole problem 6
//
// This problem is the same as problem 5, but with temperature that varies 
// in the x direction.
//
//	Problem equations:
//
//	Domaine:
//
//	Initial conditions and BC's:
//
//	Solution: 
//
//*****************************************************************************
void moleProblem6(int myid){
}

//*****************************************************************************
// Mole problem 7
//
// This problem models the production, flow, and decay of I-135 in a 1D pipe. 
// The production term is driven by an externally coupled fission source. 
// The problem should be simulated for 48 hours.
//
//	Problem equations:
//
//	Domaine:
//
//	Initial conditions and BC's:
//
//	Solution: 
//
//*****************************************************************************
void moleProblem7(int myid){
}

//*****************************************************************************
// Mole problem 8
//
// This problem builds on Problem 7 by generating Xe-135 from the I-135 decay 
// as well as from fission. Xe-135 absorption of neutrons is also modeled 
// while both I-135 and Xe-135 are transported around the loop.
//
//	Problem equations:
//
//	Domaine:
//
//	Initial conditions and BC's:
//
//	Solution: 
//
//*****************************************************************************
void moleProblem8(int myid){
}

//*****************************************************************************
// Mole problem 9
//
// This problem builds on Problem 8 by allowing Xe-135 to enter and grow 
// pre-existing gas bubbles. These bubbles are assumed to have some initial 
// distribution made entire of helium.
//
//	Problem equations:
//
//	Domaine:
//
//	Initial conditions and BC's:
//
//	Solution: 
//
//*****************************************************************************
void moleProblem9(int myid){
}

//*****************************************************************************
// Mole problem 10
//
// This problem models the primary uranium isotopes along with Pu-239 and 
// Np-239. The isotopes all start with a flat distribution and decay without 
// moving.
//
//	Problem equations:
//		dCi/dt = sum^{9}_{j=1} A_{i,j}*C_{j}
//
//		i = 1, U-233
//			 2, U-234
//			 3, U-235
//			 4, U-236
//			 5, U-237
//			 6, U-238
//			 7, U-239
//			 8, Pu-239
//			 9, Np-239
//
//	Domaine:
//		x = [0, 400]	cm
//		t = [0, 10]		y
//
//	Initial conditions and BC's:
//		Ci(x, 0) = 1e10
//
//	Solution: 
//		Computed using matlab
//
//*****************************************************************************
void moleProblem10(int myid){
	double t;
	int steps = 10;
	double depletionTime = 500.; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps;
	int xCells = 1, yCells = 1;
	double xLength = 1.0, yLength = 4.0;
	MatrixD anaSolution;
	std::vector<int> ids;
	std::string path = getDataPath();
	std::string outputFileName = "moleProblem10.out";
	std::vector<std::string> solvers {"CRAM", "hyperbolic", "parabolic"};
	//"pade-method1", "pade-method2", "taylor"};
	const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

	// File names for setting up problems
	std::string speciesNamesFile = path + "MoleP10SpeciesInputNames.dat";
	std::string speciesDecayFile = path + "MoleP10SpeciesInputDecay.dat";
	std::string solutionfname = path + "moleP10Solution.csv";

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	// Builds species object
	speciesDriver spec = speciesDriver(&model);
	// Sets the neutron flux
	model.setSystemNeutronFlux(1.e13);
	// Reads in the analytical solution
	readCSV(anaSolution, solutionfname);

	// Loop over matrix exp solvers
	for (std::string &solverType : solvers){
		std::ofstream outputFile;
		// Adds the species 
		ids = spec.addSpeciesFromFile(speciesNamesFile);
		// Set source terms
		spec.setSpeciesSourceFromFile(speciesDecayFile);
		// Write transition matrix to file
		spec.writeTransitionMatrixToFile("transitionMatrixMoleP10.csv");
		// Opens outputfile
		if (myid == 0){
			outputFile.open(outputFileName, std::ios_base::app);
			outputFile.precision(16); 
			// Write solver name
			outputFile << "solverName: " << solverType << std::endl;
		}
		// name of the csv solution output file
		std::string outputFileNameCSV = "moleProblem10"+solverType+"Steps10.csv";
		// Sets the matrix exp solver
		spec.setMatrixExpSolver(solverType);
		// sets init time
		t = 0.0;
		// sets the soltuion data array
		MatrixD solData = MatrixD::Zero(10,10);

		// Loop over the solve tiems
		for (int k = 0; k < steps; k++){
			t = t + dt;
			spec.solve(t);
			if (myid == 0){
			outputFile << "Time: " << t << std::endl;
			// Loop over species to print solution
				for (int id = 0; id < ids.size(); id++){
					std::string name = spec.getSpeciesName(0, 0, ids[id]);
					double con = spec.getSpecies(0, 0, ids[id]);
					outputFile << name << " " << con << std::endl;
					solData(id, k) = con;
				}
			}
		}
		if (myid == 0){
			outputFile << " " << std::endl;
			// Write out solution to matrix 
			writeCSV(solData, outputFileNameCSV);
			double error = computeRelativeRMSE(anaSolution, solData);
			std::cout << solverType << " " << error << std::endl;
			assert(error < 1.e-12);
		}
		// Cleans species
		spec.clean();
	}
}

//*****************************************************************************
// Mole problem 11
//
// This problem builds on Problem 10 by adding transport of the isotopes 
// around the loop.
//
//	Problem equations:
//
//	Domaine:
//
//	Initial conditions and BC's:
//
//	Solution: 
//
//*****************************************************************************
void moleProblem11(int myid){
}

//*****************************************************************************
// Mole problem 12
//
// This problem is the same as Problem 10, but includes additional transitions 
// due to a coupled neutron flux. The transition matrix was built using the 
// neutron flux shown below, in pyLibowski.
//
//	Problem equations:
//		dCi/dt = sum^{9}_{j=1} A_{i,j}*C_{j}
//
//		i = 1, U-233
//			 2, U-234
//			 3, U-235
//			 4, U-236
//			 5, U-237
//			 6, U-238
//			 7, U-239
//			 8, Pu-239
//			 9, Np-239
//
//	Domaine:
//		x = [0, 400]	cm
//		t = [0, 10]		y
//		phi = 1e13		1/cm^2/s
//
//	Initial conditions and BC's:
//		Ci(x, 0) = 1e10
//
//	Solution: 
//		Computed using matlab
//
//*****************************************************************************
void moleProblem12(int myid){
	double t;
	int steps = 10;
	double depletionTime = 500.; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps;
	int xCells = 1, yCells = 1;
	double xLength = 1.0, yLength = 4.0;
	MatrixD anaSolution;
	std::vector<int> ids;
	std::string path = getDataPath();
	std::string outputFileName = "moleProblem12.out";
	std::vector<std::string> solvers {"CRAM", "hyperbolic", "parabolic"};
	//"pade-method1", "pade-method2", "taylor"};

	// File names for setting up problems
	std::string speciesNamesFile = path + "MoleP12SpeciesInputNames.dat";
	std::string speciesDecayFile = path + "MoleP12SpeciesInputDecay.dat";
	std::string speciesTransFile = path + "MoleP12SpeciesInputTrans.dat";
	std::string solutionfname = path + "moleP12Solution.csv";

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	// Builds species object
	speciesDriver spec = speciesDriver(&model);
	// Sets the neutron flux
	model.setSystemNeutronFlux(1.e13);
	// Reads in the analytical solution
	readCSV(anaSolution, solutionfname);

	// Loop over matrix exp solvers
	for (std::string &solverType : solvers){
		std::ofstream outputFile;
		// Adds the species 
		ids = spec.addSpeciesFromFile(speciesNamesFile);
		// Set source terms
		spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);
		// Write transition matrix to file
		spec.writeTransitionMatrixToFile("transitionMatrixMoleP12.csv");
		// Opens outputfile
		if (myid == 0){
			outputFile.open(outputFileName, std::ios_base::app);
			outputFile.precision(16); 
			// Write solver name
			outputFile << "solverName: " << solverType << std::endl;
		}
		// name of the csv solution output file
		std::string outputFileNameCSV = "moleProblem12"+solverType+"Steps10.csv";
		// Sets the matrix exp solver
		spec.setMatrixExpSolver(solverType);
		// sets init time
		t = 0.0;
		// sets the soltuion data array
		MatrixD solData = MatrixD::Zero(10,10);

		// Loop over the solve tiems
		for (int k = 0; k < steps; k++){
			t = t + dt;
			spec.solve(t);
			if (myid == 0){
				outputFile << "Time: " << t << std::endl;
				// Loop over species to print solution
				for (int id = 0; id < ids.size(); id++){
					std::string name = spec.getSpeciesName(0, 0, ids[id]);
					double con = spec.getSpecies(0, 0, ids[id]);
					outputFile << name << " " << con << std::endl;
					solData(id, k) = con;
				}
			}
		}
		if (myid == 0){
			outputFile << " " << std::endl;
			// Write out solution to matrix 
			writeCSV(solData, outputFileNameCSV);
			double error = computeRelativeRMSE(anaSolution, solData);
			std::cout << solverType << " " << error << std::endl;
			assert(error < 1.e-12);
		}
		// Cleans species
		spec.clean();
	}
}

//*****************************************************************************
// Mole problem 13
//
// This problem is the same as Problem 12, but includes transport due to 
// velocity of the fluid.
//
//	Problem equations:
//
//	Domaine:
//
//	Initial conditions and BC's:
//
//	Solution: 
//
//*****************************************************************************
void moleProblem13(int myid){
}

int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	moleProblem1(myid); 
	moleProblem2(myid);
	moleProblem3(myid);
	
	//moleProblem4(myid);
	//moleProblem5(myid);
	//moleProblem6(myid);
	//moleProblem7(myid);
	//moleProblem8(myid);
	//moleProblem9(myid);
	moleProblem10(myid);
	//moleProblem11(myid);
	moleProblem12(myid);
	//moleProblem13(myid);

	mpi.finalize();
}
