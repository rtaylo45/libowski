#include <string>
#include <math.h>
#include <cmath>
#include <iostream>

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
//			lambda = 1.0
//
//	Initial conditions and BC's:
//			C(x, 0) = x
//			C(0,t) = C(100,t)
//
//	Solution:
//			C(x,t) = x*e^(-lambda*t)
//
//*****************************************************************************
void moleProblem1(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10, 20};
	std::vector<double> steps = {1, 2, 5, 10};
	double xLength = 100, yLength = 0.0; // cm
	double tEnd = 20.0;	// seconds
	double lambda = 1.0;	// 1/s
	double cSol, cCon;
	int cID;
	double x1, x2, initCon, dx, xc, dt, t;
	double linfError = 0.0;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> cCoeffs = {-lambda};
	// Sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "moleProblem1.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "dx"	<< " " << "dt" << " " << "Max l2 Error" << "\n";

	// Loops over number of cells
	for (int &xCells : numOfxCells){
		// Build the Mesh
		modelMesh model(xCells, yCells, xLength, yLength);
		// Build species driver
		speciesDriver spec = speciesDriver(&model);


		// Loops over number of time steps
		for (double &numOfSteps	: steps){
			dt = tEnd/numOfSteps;

			// Add specs
			cID = spec.addSpecies(1.0, 0.0, 0.0);

			// Sets the intial condition and sources
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	

					// Calculates the x positions as the cell faces
					dx = cell->dx;
					xc = cell->x;
					x2 = xc + dx/2;
					x1 = xc - dx/2;

					// Calculates the initial concentration from MVT. 
					initCon = (1./(2.*dx))*(x2*x2 - x1*x1);		

					spec.setSpeciesCon(i, j, cID, initCon);

					// Sets the sources
					spec.setSpeciesSource(i, j, cID, cCoeffs, 0.0);
				}
			}
			t = 0.0;
			// Solve the problem
			for (int step = 1; step <= numOfSteps; step++){
				t = step*dt;
				// Solve with CRAM
				spec.solve(t);

				// Gets species Concentrations
				if (myid==0){
					for (int i = 0; i < xCells; i++){
						for (int j = 0; j < yCells; j++){
							cell = model.getCellByLoc(i,j);	

							// Caclulate analytical solution
							xc = cell->x;
							cSol = xc*exp(-lambda*t);

							// Get libowski solution
							cCon = spec.getSpecies(i, j, cID);

							//assert(isApprox(VSol, VCon, 1e-5, 1e-4));
							linfError = std::max(linfError, std::abs(cSol-cCon)/cSol);
						}
					}
				}
				outputFile << dx << " " << dt << " " << linfError << "\n";
			}
			// Clean species
			spec.clean();
		}
		spec.clean();
		model.clean();
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
//			lambda = 1.0
//
//	Initial conditions and BC's:
//			C(x, 0) = x
//			C(0,t) = C(100,t)
//
//	Solution:
//			C(x,t) = (x-v*t)*e^(-lambda*t)
//
//*****************************************************************************
void moleProblem2(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{5};
	std::vector<double> steps = {20};
	double xLength = 100, yLength = 0.0; // cm
	double tEnd = 20.0;	// seconds
	double lambda = 0.1;	// 1/s
	double velocity = 2.0; // cm/s
	double cSol, cCon;
	int cID;
	double x1, x2, initCon, dx, xc, dt, t;
	double linfError = 0.0;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> ccoeffs = {-lambda};
	// sets the ouput file name
	std::ofstream outputFile;
	outputFileName = "moleproblem2.out";
	outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
	outputFile << "dx"	<< " " << "dt" << " " << "max l2 error" << "\n";

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

			// add specs
			cID = spec.addSpecies(1.0, 0.0, 0.0);

			// set periodic bcs
			spec.setBoundaryCondition("periodic","east", cID);
			spec.setBoundaryCondition("periodic","west", cID);

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
					initCon = (1./(2.*dx))*(x2*x2 - x1*x1);		

					spec.setSpeciesCon(i, j, cID, initCon);

					// sets the sources
					spec.setSpeciesSource(i, j, cID, ccoeffs, 0.0);
				}
			}
			t = 0.0;
			// solve the problem
			for (int step = 1; step <= numofsteps; step++){
				t = step*dt;
				// solve with cram
				spec.solve(t);
				//spec.solveimplicit(t);

				std::cout << t << std::endl;
				// gets species concentrations
				if (myid==0){
					for (int i = 0; i < xCells; i++){
						for (int j = 0; j < yCells; j++){
							cell = model.getCellByLoc(i,j);	

							// caclulate analytical solution
							xc = cell->x;
							cSol = (xc - velocity*t)*exp(-lambda*t);

							// get libowski solution
							cCon = spec.getSpecies(i, j, cID);

							//assert(isapprox(vsol, vcon, 1e-5, 1e-4));
							linfError = std::max(linfError, std::abs(cSol-cCon)/cSol);
							std::cout << cSol << " " << cCon << std::endl;
						}
					}
				}
				//outputfile << dx << " " << dt << " " << linferror << "\n";
				std::cout << dx << " " << dt << " " << linfError << "\n";
			}
			// clean species
			spec.clean();
		}
		spec.clean();
		model.clean();
	}
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
//			lambda = 1.0
//			R(T) = A1*e^(A2*T)
//			A1 = 1.0
//			A2 = 1./3.
//			T = 900
//
//	Initial conditions and BC's:
//			Ci(x, 0) = x
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
	double A1 = 1.0, A2 = 5000./8.31, T = 900.;
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
	std::cout << "R " << R << std::endl;

	// loops over number of cells
	for (int &xCells : numOfxCells){
		// build the mesh
		modelMesh model(xCells, yCells, xLength, yLength);
		// build species driver
		speciesDriver spec = speciesDriver(&model);

		// loops over number of time steps
		for (double &numofsteps	: steps){
			dt = tEnd/numofsteps;

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
					initCon = (1./(2.*dx))*(x2*x2 - x1*x1);		

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

				std::cout << t << std::endl;
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
							std::cout << CiSol << " " << CiCon << " " << CwSol << " " << CwCon << std::endl;
						}
					}
				}
				//outputfile << dx << " " << dt << " " << linferror << "\n";
				std::cout << dx << " " << dt << " " << linfError << "\n";
			}
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
//			lambda = 1.0
//			v = 2.0
//			R(T) = A1*e^(A2*T)
//			A1 = 1.0
//			A2 = 1./3.
//			T = 900
//
//	Initial conditions and BC's:
//			Ci(x, 0) = x
//			Ci(0,t) = Ci(100,t)
//			Cw(x, 0) = 0
//			Cw(0,t) = Cw(100,t)
//
//	Solution: Need to fix this
//			Ci(x,t) = Ci0 - R(T)*t
//			Cw(x,t) = R(T)*t
//
//*****************************************************************************
void moleProblem4(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{10};
	std::vector<double> steps = {1};
	double xLength = 100, yLength = 0.0; // cm
	double tEnd = 20.0;	// seconds
	double A1 = 1.0, A2 = 5000./8.31, T = 900.;
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
	outputFile << "dx"	<< " " << "dt" << " " << "max l2 error" << "\n";

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
					initCon = (1./(2.*dx))*(x2*x2 - x1*x1);		

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

				std::cout << t << std::endl;
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
							std::cout << CiSol << " " << CiCon << " " << CwSol << " " << CwCon << std::endl;
						}
					}
				}
				//outputfile << dx << " " << dt << " " << linferror << "\n";
				std::cout << dx << " " << dt << " " << linfError << "\n";
			}
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

	//moleProblem1(myid);
	//moleProblem2(myid);
	//moleProblem3(myid);
	moleProblem4(myid);


}
