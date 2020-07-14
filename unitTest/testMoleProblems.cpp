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
//			x = [0, 100]
//			t = [0, 20]
//			lambda = 0.01
//
//	Initial conditions and BC's:
//			C(x, 0) = 1000.(x + 1)
//			C(0,t) = C(100,t)
//
//	Solution:
//			C(x,t) = C(x,0)*e^(-lambda*t)
//
//*****************************************************************************
void moleProblem1(int myid){
	int yCells = 1, xCells = 1000;
	std::vector<double> steps = {1, 2, 4, 8, 20, 40, 80, 200, 400};
	std::vector<std::string> solvers {"hyperbolic", "pade-method2", "taylor"};
	double xLength = 100., yLength = 0.0; // cm
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
						outputFile << xc << " " << cCon << "\n";
						assert(isApprox(cCon, cSol));
					}
				}
				outputFile << "\n";
				printf("%15s %4.2f %4.2e %4.2e %3.5f \n", solverType.c_str(), dt, relativeError,
						std::pow(rmse/float(xCells), 0.5), duration.count()/1.e6);
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
//			x = [0, 100]
//			t = [0, 20]
//			v = 2.0
//			lambda = 0.01
//
//	Initial conditions and BC's:
//			C(x, 0) = 1000.0
//			C(0,t) = 1000.0
//
//	Solution:
//			C(x,t) = TBD
//
//*****************************************************************************
void moleProblem2(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10, 100, 1000};
	//std::vector<double> steps = {5};
	std::vector<double> steps = {1, 2, 4, 8, 20, 40, 80, 200, 400};
	//std::vector<double> steps = {1, 2, 4, 8, 20, 40};
	std::vector<std::string> solvers {"hyperbolic","pade-method2", "taylor"};
	//std::vector<std::string> solvers {"hyperbolic"};
	double xLength = 100, yLength = 0.0; // cm
	double tEnd = 20.0;	// seconds
	double lambda = 0.01;	// 1/s
	//double lambda = 0.0;	// 1/s
	double velocity = 2.0; // cm/s
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
				spec.setBoundaryCondition("free flow","east", cID);

				// sets the intial condition and sources
				for (int i = 0; i < xCells; i++){
					for (int j = 0; j < yCells; j++){
						cell = model.getCellByLoc(i,j);	

						// calculates the x positions as the cell faces
						dx = cell->dx;
						xc = cell->x;
						x2 = xc + dx/2.;
						x1 = xc - dx/2.;

						// calculates the initial concentration from mvt. 
						initCon = 1000.;
						//initCon = 0.0;

						spec.setSpeciesCon(i, j, cID, initCon);

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
							//std::cout << xc << std::endl; // << cSol << " " << cCon << std::endl;
							relativeError = std::abs(cSol-cCon)/cSol;
							maxRelativeError = std::max(maxRelativeError, relativeError);
							rmse += std::pow(relativeError,2.);
						}
					}
					outputFile << "\n";
					//std::cout << solverType << " " << dx << " " << dt 
					//	<< " " << percentError << "\n";
					// clean species
					printf("%15s %2.3f %4.2E %4.2e %4.2e %3.5f \n", solverType.c_str(), dt, dx, 
							maxRelativeError, std::pow(rmse/float(xCells), 0.5), duration.count()/1.e6);
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
//			dCi/dt = -R(T)
//			dCw/dt = R(T)
//
//	Domaine:
//			x = [0, 100]
//			t = [0, 20]
//			R(T) = A1*e^(A2/T)
//			A1 = 0.1
//			A2 = -500.
//			T = 900
//
//	Initial conditions and BC's:
//			Ci(x, 0) = 1000.0
//			Ci(0,t) = Ci(100,t)
//			Cw(x, 0) = 0
//			Cw(0,t) = Cw(100,t)
//
//	Solution:
//			Ci(x,t) = Ci0 - R(T)*t
//			Cw(x,t) = R(T)*t
//
//*****************************************************************************
void moleProblem3(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10};
	std::vector<double> steps = {1};
	double xLength = 100, yLength = 0.0; // cm
	double tEnd = 20.0;	// seconds
	double A1 = 1.0, A2 = 500, T = 900.;
	double R = A1*exp(-A2/T);
	double CiSol, CwSol, CiCon, CwCon;
	int CiID, CwID;
	double x1, x2, initCon, dx, xc, dt, t;
	double linfError = 0.0;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> Cicoeffs = {0.0, 0.0};
	std::vector<double> Cwcoeffs = {0.0, 0.0};
	// sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "moleproblem3.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "dx"	<< " " << "dt" << " " << "max l2 error" << "\n";
	outputFile << "Total problem time: " << tEnd << "\n";

	// loops over number of cells
	for (int &xCells : numOfxCells){
		// build the mesh
		modelMesh model(xCells, yCells, xLength, yLength);
		// build species driver
		speciesDriver spec = speciesDriver(&model);

		// loops over number of time steps
		for (double &numofsteps	: steps){
			dt = tEnd/numofsteps;
			outputFile << "dx: " << xLength/(double)xCells << "\n";
			outputFile << "dt: " << dt << "\n";
			outputFile << "x " << "Ci " << "Cw "<< "\n";

			// add specs
			CiID = spec.addSpecies(1.0, 0.0, 0.0);
			CwID = spec.addSpecies(1.0, 0.0, 0.0);

			// sets the intial condition and sources
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	

					// calculates the x positions as the cell faces
					dx = cell->dx;
					xc = cell->x;
					x2 = xc + dx/2;
					x1 = xc - dx/2;

					// calculates the initial concentration from mvt. 
					initCon = 1000.0;

					spec.setSpeciesCon(i, j, CiID, initCon);

					// sets the sources
					spec.setSpeciesSource(i, j, CiID, Cicoeffs, -R);
					spec.setSpeciesSource(i, j, CwID, Cwcoeffs, R);
				}
			}
			t = 0.0;
			// solve the problem
			for (int step = 1; step <= numofsteps; step++){
				t = step*dt;
				// solve with cram
				spec.solve(t);
			}
			// gets species concentrations
			if (myid==0){
				for (int i = 0; i < xCells; i++){
					for (int j = 0; j < yCells; j++){
						cell = model.getCellByLoc(i,j);	

						// caclulate analytical solution
						xc = cell->x;
						CiSol = xc - R*t + 1000.;
						CwSol = R*t;

						// get libowski solution
						CiCon = spec.getSpecies(i, j, CiID);
						CwCon = spec.getSpecies(i, j, CwID);

						//assert(isapprox(vsol, vcon, 1e-5, 1e-4));
						//linfError = std::max(linfError, std::abs(csol-ccon)/csol);
						outputFile << xc << " " << CiCon << " " << CwCon << "\n";
					}
				}
			}
			outputFile << "\n";
			//std::cout << dx << " " << dt << " " << linfError << "\n";
			// clean species
			spec.clean();
		}
		spec.clean();
		model.clean();
	}
}

//*****************************************************************************
// Mole problem 4
//
// This problem is the same as problem 3, but with the liquid flowing at a 
// constant velocity in the positive x direction
//
//	Problem equations:
//			dCi/dt = -v*dCi/dx-R(T)
//			dCw/dt = R(T)
//
//	Domaine:
//			x = [0, 100]
//			t = [0, 20]
//			v = 2.0
//			R(T) = A1*e^(A2/T)
//			A1 = 0.1
//			A2 = -500.
//			T = 900
//
//	Initial conditions and BC's:
//			Ci(x, 0) = 1000.0
//			Ci(0,t) = 1000.0
//			Cw(x, 0) = 0
//			Cw(0,t) = Cw(100,t)
//
//	Solution: 
//			Ci(x,t) = TBD
//			Cw(x,t) = TBD
//
//*****************************************************************************
void moleProblem4(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10};
	std::vector<double> steps = {1};
	double xLength = 100, yLength = 0.0; // cm
	double tEnd = 10.0;	// seconds
	double A1 = 0.1, A2 = 500., T = 900.;
	double velocity = 2.0; // cm/s
	double R = A1*exp(-A2/T);
	double CiSol, CwSol, CiCon, CwCon;
	int CiID, CwID;
	double x1, x2, initCon, dx, xc, dt, t;
	double linfError = 0.0;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> Cicoeffs = {0.0, 0.0};
	std::vector<double> Cwcoeffs = {0.0, 0.0};
	// sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "moleproblem4.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "Total problem time: " << tEnd << "\n";

	// loops over number of cells
	for (int &xCells : numOfxCells){
		// build the mesh
		modelMesh model(xCells, yCells, xLength, yLength);
		// build species driver
		speciesDriver spec = speciesDriver(&model);
		// set x velocity
		model.setConstantXVelocity(velocity);

		// loops over number of time steps
		for (double &numofsteps	: steps){
			dt = tEnd/numofsteps;
			outputFile << "dx: " << xLength/(double)xCells << "\n";
			outputFile << "dt: " << dt << "\n";
			outputFile << "x " << "Ci " << "Cw "<< "\n";

			// add specs
			CiID = spec.addSpecies(1.0, 0.0, 0.0);
			CwID = spec.addSpecies(1.0, 0.0, 0.0);

			// sets the intial condition and sources
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	

					// calculates the x positions as the cell faces
					dx = cell->dx;
					xc = cell->x;
					x2 = xc + dx/2;
					x1 = xc - dx/2;

					// calculates the initial concentration from mvt. 
					initCon = (5.*(std::pow(M_PI, 0.5))/(dx))*(erf(3.-x1/10.) - erf(3.-x2/10.));
					initCon += (1/dx)*(10.*x2 - 10.*x1);

					spec.setSpeciesCon(i, j, CiID, initCon);

					// sets the sources
					spec.setSpeciesSource(i, j, CiID, Cicoeffs, -R);
					spec.setSpeciesSource(i, j, CwID, Cwcoeffs, R);
				}
			}
			t = 0.0;
			// solve the problem
			for (int step = 1; step <= numofsteps; step++){
				t = step*dt;
				// solve with cram
				spec.solve(t);
			}
			// gets species concentrations
			if (myid==0){
				for (int i = 0; i < xCells; i++){
					for (int j = 0; j < yCells; j++){
						cell = model.getCellByLoc(i,j);	

						// caclulate analytical solution
						xc = cell->x;
						CiSol = 0.0;
						CwSol = R*t;

						// get libowski solution
						CiCon = spec.getSpecies(i, j, CiID);
						CwCon = spec.getSpecies(i, j, CwID);

						//assert(isapprox(vsol, vcon, 1e-5, 1e-4));
						//linfError = std::max(linfError, std::abs(csol-ccon)/csol);
						outputFile << xc << " " << CiCon << " " << CwCon << "\n";
					}
				}
			}
			outputFile << "\n";
			//std::cout << dx << " " << dt << " " << linfError << "\n";
			// clean species
			spec.clean();
		}
		spec.clean();
		model.clean();
	}
}

//*****************************************************************************
// Mole problem 5
//
// This problem is the same as problem 4 but with temperature that varies in 
// the x direction. Source term is constant but varies in the x direction,
// so need to use MVT and integrate over the cell. 
//
//	Problem equations:
//			dCi/dt = -v*dCi/dx-R(T)
//			dCw/dt = R(T)
//
//	Domaine:
//			x = [0, 100]
//			t = [0, 20]
//			v = 2.0
//			R(T) = A1*e^(A2/T)
//			T = 850 + x
//			A1 = 0.1
//			A2 = -500.
//
//	Initial conditions and BC's:
//			Ci(x, 0) = 1000.0
//			Ci(0,t) = 1000.0
//			Cw(x, 0) = 0
//			Cw(0,t) = Cw(100,t)
//
//	Solution: 
//			C(x,t) = TBD
//			Cw(x,t) = TBD
//
//*****************************************************************************
void moleProblem5(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10};
	std::vector<double> steps = {1};
	double xLength = 100, yLength = 0.0; // cm
	double tEnd = 10.0;	// seconds
	double A1 = 0.1, A2 = 500., T = 900.;
	double velocity = 2.0; // cm/s
	double R = A1*exp(-A2/T);
	double CiSol, CwSol, CiCon, CwCon;
	int CiID, CwID;
	double x1, x2, initCon, dx, xc, dt, t;
	double linfError = 0.0;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> Cicoeffs = {0.0, 0.0};
	std::vector<double> Cwcoeffs = {0.0, 0.0};
	// sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "moleproblem5.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "Total problem time: " << tEnd << "\n";

	// loops over number of cells
	for (int &xCells : numOfxCells){
		// build the mesh
		modelMesh model(xCells, yCells, xLength, yLength);
		// Add BC surface
		model.addBoundarySurface("west");
		// build species driver
		speciesDriver spec = speciesDriver(&model);
		// set x velocity
		model.setConstantXVelocity(velocity);

		// loops over number of time steps
		for (double &numofsteps	: steps){
			dt = tEnd/numofsteps;
			outputFile << "dx: " << xLength/(double)xCells << "\n";
			outputFile << "dt: " << dt << "\n";
			outputFile << "x " << "Ci " << "Cw "<< "\n";

			// add specs
			CiID = spec.addSpecies(1.0, 0.0, 0.0);
			CwID = spec.addSpecies(1.0, 0.0, 0.0);

			// sets the intial condition and sources
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	

					// calculates the x positions as the cell faces
					dx = cell->dx;
					xc = cell->x;
					x2 = xc + dx/2;
					x1 = xc - dx/2;

					// calculates the initial concentration from mvt. 
					initCon = (5.*(std::pow(M_PI, 0.5))/(dx))*(erf(3.-x1/10.) - erf(3.-x2/10.));
					initCon += (1/dx)*(10.*x2 - 10.*x1);

					spec.setSpeciesCon(i, j, CiID, initCon);

					// sets the sources
					spec.setSpeciesSource(i, j, CiID, Cicoeffs, -R);
					spec.setSpeciesSource(i, j, CwID, Cwcoeffs, R);
				}
			}
			t = 0.0;
			// solve the problem
			for (int step = 1; step <= numofsteps; step++){
				t = step*dt;
				// solve with cram
				spec.solve(t);
			}
			// gets species concentrations
			if (myid==0){
				for (int i = 0; i < xCells; i++){
					for (int j = 0; j < yCells; j++){
						cell = model.getCellByLoc(i,j);	

						// caclulate analytical solution
						xc = cell->x;
						CiSol = xc - R*t;
						CwSol = R*t;

						// get libowski solution
						CiCon = spec.getSpecies(i, j, CiID);
						CwCon = spec.getSpecies(i, j, CwID);

						//assert(isapprox(vsol, vcon, 1e-5, 1e-4));
						//linfError = std::max(linfError, std::abs(csol-ccon)/csol);
						outputFile << xc << " " << CiCon << " " << CwCon << "\n";
					}
				}
			}
			outputFile << "\n";
			//std::cout << dx << " " << dt << " " << linfError << "\n";
			// clean species
			spec.clean();
		}
		spec.clean();
		model.clean();
	}
}
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	moleProblem1(myid); moleProblem2(myid);
	//moleProblem3(myid);
	//moleProblem4(myid);
	//moleProblem5(myid);

	mpi.finalize();
}
