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
// Problem problem 1
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
// Problem problem 2
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
//
//	Initial conditions and BC's:
//			C(x, 0) = x
//			C(0,t) = C(100,t)
//
//	Solution:
//			C(x,t) = x*e^(-lambda*t)
//
//*****************************************************************************
void moleProblem2(int myid){
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

int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	moleProblem1(myid);
	//moleProblem2(myid);


}
