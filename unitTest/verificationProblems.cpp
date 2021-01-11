#define _USE_MATH_DEFINES
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>
#include <unsupported/Eigen/CXX11/Tensor>
#include <assert.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iomanip>

#include "mpiProcess.h"
#include "modelMesh.h"
#include "speciesDriver.h"
#include "meshCell.h"
#include "species.h"
#include "utilBase.h"
#include "matrixTypes.h"

//*****************************************************************************
// Precursor analytical solution
//*****************************************************************************
double precursorAnalitical(double y, double vel, double a, double length, 
	double lambda){
	double sol = 0.0;
	double bottom1, bottom2;

	//bottom1 = pow(length*lambda,2.0) + pow(M_PI*vel,2.0);
	//bottom2 = pow(lambda/vel,2.0) + pow(M_PI/length,2.0);

	//sol = length*a*vel*M_PI*exp(-lambda*y/vel)/bottom1 - 
	//	(a/vel)*((M_PI/length)*cos(M_PI*y/length) - 
	//	(lambda/vel)*sin(M_PI*y/length))/bottom2;

	if (vel != 0.0){
		// if vel is not zero then y is the direction location
		sol = (a - a*exp(-lambda*y/vel))/lambda;
	}
	else{
		// if vel is zero then y is t for time
		sol = (a - a*exp(-lambda*y))/lambda;
	}

	return sol;
}

//*****************************************************************************
// Three speceis test problem
//
// dN1/dt = -lambda1*N1 + lambda3*N3
// dN2/dt = -lambda2*N2 + lambda1*N1
// dN3/dt = -lambda3*N3 + lambda2*N2
//
//	Note: The units for this problem kinda dont matter either. This is just an
//	internal ODE test, the ODE solution is provided and calculated using the 
//	same problem units. The solution is check by a python script
//*****************************************************************************
void testProblem1(int myid){
	int xCells = 1, yCells = 1;
	double xLength = 1.0, yLength = 1.0;
	double N1InitCon = 10000.0, N2InitCon = 0.0, N3InitCon = 0.0;
	double N1MM = 1.0, N2MM = 1.0, N3MM = 1.0;
	std::vector<double> steps = {600, 500, 400, 350, 300,
		250, 200, 150, 100, 75, 50, 25, 20, 15, 10, 5, 4, 3, 2, 1};
	std::vector<std::string> solvers {"CRAM", "parabolic", "hyperbolic",
		"pade-method1", "pade-method2", "taylor"};
	//std::vector<std::string> solvers {"BDF1", "BDF2", "BDF3", "BDF4", "BDF5", "BDF6"};
	double tEnd = 600.0;
	double t, dt;
	double lambda1 = 1.0/1.0e2, lambda2 = 0.5/1.0e2, lambda3 = 3.0/1.0e2;
   //double lambda1 = 1.0/1.0e2, lambda2 = 0.5/1.0e2, lambda3 = lambda1+lambda2;
	double DN1 = 0.0, DN2 = 0.0, DN3 = 0.0;
	int N1ID, N2ID, N3ID;
	double N1Con, N2Con, N3Con;
	std::vector<double> N1Coeffs = {-lambda1, 0.0, lambda3};
	std::vector<double> N2Coeffs = {lambda1, -lambda2, 0.0};
	std::vector<double> N3Coeffs = {0.0, lambda2, -lambda3};
	std::string outputFileName;

	modelMesh model(xCells, yCells, xLength, yLength);
	speciesDriver spec = speciesDriver(&model);

	// Loops over different solvers
	for (std::string &solverType : solvers){
		// Sets the species matrix exp solver
		spec.setMatrixExpSolver(solverType);
		// Sets the ouput file name
		std::ofstream outputFile;
		outputFileName = "problem1"+solverType+".out";
		outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
		outputFile << "TotalTime: "+std::to_string(tEnd)+"\n";

		// Loops over number of time steps
		for (double &numOfSteps	: steps){
			dt = tEnd/numOfSteps;
			outputFile << "TimeStepSize: "+std::to_string(dt)+"\n";

			// Sets the problem
			N1ID = spec.addSpecies(N1MM, N1InitCon, DN1);
			N2ID = spec.addSpecies(N2MM, N2InitCon, DN2);
			N3ID = spec.addSpecies(N3MM, N3InitCon, DN3);

			// Set source
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					spec.setSpeciesSource(i, j, N1ID, N1Coeffs);
					spec.setSpeciesSource(i, j, N2ID, N2Coeffs);
					spec.setSpeciesSource(i, j, N3ID, N3Coeffs);
				}
			}
			
			t = 0.0;
			// Solves the problem for each time step
			for (int step = 1; step <= numOfSteps; step++){
				t = step*dt;
				// Solve with CRAM
				spec.solve(t);

				// Gets species Concentrations
				if (myid==0){
					for (int i = 0; i < xCells; i++){
						for (int j = 0; j < yCells; j++){
							N1Con = spec.getSpecies(i, j, N1ID);
							N2Con = spec.getSpecies(i, j, N2ID);
							N3Con = spec.getSpecies(i, j, N3ID);
							//std::cout << t << " " << N1Con //<< std::endl;
							outputFile <<  std::setprecision(16) << t << " " << N1Con 
							<< " " << N2Con << " " << N3Con <<  std::endl;
						}
					}
				}
			}
			// Clean species
			spec.clean();
		}
	}
	
	model.clean();
	spec.clean();
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
void testProblem2(int myid){
	int yCells = 1;
	std::vector<int> numOfxCells{100, 150, 200};
	//std::vector<int> numOfxCells{10, 20, 30, 40, 50, 100,
	//	150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 
	//	650, 700, 750, 800, 850, 900, 950, 1000};
	double xLength = M_PI/2., yLength = 0.0;
	double numOfSteps = 1;
	double tEnd = 1.0;
	double dt = tEnd/numOfSteps, t = 0;
	double a = 0.1, b = 0.01, d = 1.0;	// Problem 2a
	//double a = 2.0, b = 1.0, d = 0.001;		// Problem 2b
	//double a = 100.0, b = 1.0, d = 0.001;		// Problem 2c
	double UCon, VCon, USol, VSol;
	int UID, VID;
	double x1, x2, xc, initCon, x, dx;
	double absError, linfError;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> Ucoeffs = {-a, 1.0};
	std::vector<double> Vcoeffs = {0.0, -b};
	std::vector<std::string> solvers {"CRAM", "parabolic", "hyperbolic", 
		"pade-method1", "pade-method2", "taylor"};

	// Loops over different solvers
	for (std::string &solverType : solvers){

		std::ofstream outputFile;
		outputFileName = "problem2"+solverType+".out";
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

			// Solve with CRAM
			auto start = std::chrono::high_resolution_clock::now();
			for (int step = 1; step <= numOfSteps; step++){
				t = step*dt;
				spec.solve(t);
			}
			auto end = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
				end - start);

			linfError = 0.0;
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
					
						assert(isApprox(USol, UCon, 1e-5, 1e-4));
						assert(isApprox(VSol, VCon, 1e-5, 1e-4));
						absError = std::abs(USol-UCon) + std::abs(VSol-VCon);
						linfError = std::max(linfError, absError);

					}
				}
			}
			outputFile << std::setprecision(16) << " " << xCells << " " << 
				dx << " " << linfError << " " 
				<< duration.count()/1.e6 << std::endl;
			std::cout << solverType << " " << xCells << " " << linfError <<
				" " << duration.count()/1.e6 << std::endl;
			assert(linfError < 1.e-5);

			model.clean();
			spec.clean();
		}
	}

}

//*****************************************************************************
// 2 species 1D diffusion problem using Krylov subspace method on the Pade
// solvers. Taken from the following paper 
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
void testProblem2Krylov(int myid){
	int yCells = 1;
	int xCells = 800;
	//std::vector<int> krylovDims = lineSpace(5, 1000, 996);
	std::vector<int> krylovDims = lineSpace(5, 100, 96);
	//std::vector<int> krylovDims = {100};
	double xLength = M_PI/2., yLength = 0.0;
	double numOfSteps = 1;
	double tEnd = 1.0;
	double dt = tEnd/numOfSteps, t = 0;
	double a = 0.1, b = 0.01, d = 1.0;	// Problem 2a
	//double a = 2.0, b = 1.0, d = 0.001;		// Problem 2b
	//double a = 100.0, b = 1.0, d = 0.001;		// Problem 2c
	double UCon, VCon, USol, VSol;
	int UID, VID;
	double x1, x2, xc, initCon, x, dx;
	double relativeErr, linfError;
	meshCell* cell = nullptr;
	//std::string outputFileName;
	std::string outputFileName = "problem2Krylov.out";
	remove(outputFileName.c_str());
	std::vector<double> Ucoeffs = {-a, 1.0};
	std::vector<double> Vcoeffs = {0.0, -b};
	std::vector<std::string> solvers {"pade-method1", "pade-method2"};
	//std::vector<std::string> solvers {"pade-method1", "pade-method2", 
	//		"taylor"};
	FILE * pOutputFile;
	pOutputFile = fopen(outputFileName.c_str(), "a");

	// Build the Mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");
	// Build species driver
	speciesDriver spec = speciesDriver(&model);

	// Loops over different solvers
	for (std::string &solverType : solvers){

		//std::ofstream outputFile;
		//outputFileName = "problem2krylov"+solverType+".out";
		//outputFile.open(outputFileName, std::ios::out | std::ios::trunc);
		fprintf(pOutputFile, "solver: %s \n", solverType.c_str());
		fprintf(pOutputFile, "elements: %ld \n", krylovDims.size());
		fprintf(pOutputFile, "%s %s %s %s \n", "variables", "dim", "linf", 
			"runtime");

		for (int &krylovDim : krylovDims){

			// Add species. I will add the initial condition later
			UID = spec.addSpecies(1.0, 0.0, d);
			VID = spec.addSpecies(1.0, 0.0, d);

			// Add BCs
			spec.setBoundaryCondition("newmann","west", UID, 0.0);
			spec.setBoundaryCondition("dirichlet","east", UID, 0.0);
			spec.setBoundaryCondition("newmann","west", VID, 0.0);
			spec.setBoundaryCondition("dirichlet","east", VID, 0.0);

			// Sets the species matrix exp solver
			spec.setMatrixExpSolver(solverType, true, krylovDim);
			//spec.setMatrixExpSolver(solverType);

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

				linfError = 0.0;
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

							if (krylovDim > 80){
								assert(isApprox(USol, UCon, 1e-6, 1e-5));
								assert(isApprox(VSol, VCon, 1e-6, 1e-5));
							}
							relativeErr = std::abs(USol-UCon) + 
								std::abs(VSol-VCon);
							linfError = std::max(linfError, relativeErr);

						}
					}
				}
				std::cout << std::setprecision(16) << solverType << " " << 
					krylovDim << " "  << linfError << " " << 
					duration.count()/1.e6 << std::endl;
				//outputFile << xCells << " " << krylovDim << " " 
				//	<< linfError << " " << duration.count()/1.e6 << std::endl;
				fprintf(pOutputFile, "%2d %8.7e %8.7e \n", krylovDim, 
					linfError, duration.count()/1.e6);

			}
		spec.clean();
		}
		fprintf(pOutputFile, "\n");
	}
	model.clean();
	fprintf(pOutputFile, "end");
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
void testProblem2IntegratorMethods(int myid){
	int yCells = 1;
	int xCells = 500;
	std::vector<double> steps = lineSpace(200.,300.,100);
	double xLength = M_PI/2., yLength = 0.0;
	double tEnd = 1.0;
	double dt, t = 0;
	double a = 0.1, b = 0.01, d = 1.0;	// Problem 2a
	//double a = 2.0, b = 1.0, d = 0.001;		// Problem 2b
	//double a = 100.0, b = 1.0, d = 0.001;		// Problem 2c
	double UCon, VCon, USol, VSol;
	int UID, VID;
	double x1, x2, xc, initCon, x, dx;
	double linfErrorU, linfErrorV;
	meshCell* cell = nullptr;
	std::string outputFileName;
	std::vector<double> Ucoeffs = {-a, 1.0};
	std::vector<double> Vcoeffs = {0.0, -b};
	std::vector<std::string> solvers {"BDF4", "BDF5", "BDF6"};

	// Build the Mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");
	// Build species driver
	speciesDriver spec = speciesDriver(&model);

	// Loops over different solvers
	for (std::string &solverType : solvers){

		std::ofstream outputFile;
		outputFileName = "problem2"+solverType+".out";
		outputFile.open(outputFileName, std::ios::out | std::ios::trunc);

		
		// Loops over number of time steps
		for (double &numOfSteps	: steps){
			dt = tEnd/numOfSteps;

			// sets the solver
			spec.setIntegratorSolver("implicit", solverType);

			// Add species. I will add the initial condition later
			UID = spec.addSpecies(1.0, 0.0, d);
			VID = spec.addSpecies(1.0, 0.0, d);

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

			auto start = std::chrono::high_resolution_clock::now();
			for (int step = 1; step <= numOfSteps; step++){
				t = step*dt;
				spec.solveImplicit(t);

			}
			auto end = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
				end - start);

			linfErrorU = 0.0;
			linfErrorV = 0.0;
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

						assert(isApprox(USol, UCon, 1e-5, 1e-4));
						assert(isApprox(VSol, VCon, 1e-5, 1e-4));
						linfErrorU = std::max(linfErrorU, std::abs(USol-UCon));
						linfErrorV = std::max(linfErrorV, std::abs(VSol-VCon));
					}
				}
			}
			//outputFile << std::setprecision(16) << " " << dt << " " << linfErrorU << " " 
			std::cout << " " << dx << " " << linfErrorU << " " 
				<< linfErrorV << " " << duration.count()/1.e6 << std::endl;

			spec.clean();
		}
		spec.clean();
	}
}

//*****************************************************************************
// Single species decay, 1D 
//
// Diff eqs:
//		dC_k/dt = -v*dC_k/dx - lambda*C_k + source
//
//		k = 1, 2, 3, 4, 5, 6
//
//	Boundary Conditions:
//		Periodic
//
//	Initial Conditions:
//		C = 0
//
//	Solution:
//
//	Note: The source terms were taken from a problem that Aaron Grahm built.
//	They should be in g/cm^3/s for the constant source terms. Need to convert 
//	this value.
//*****************************************************************************
void testProblem3(int myid){
	int xCells = 1, yCells = 50;
	double xLength = 0.0, yLength = 1.0; // m
	double v = 0.;
	//std::vector<double> steps = lineSpace(1.,300.,300);
	std::vector<double> steps = {1};
	double tEnd = 1.0;
	double t;
	double dt;
	double cCon1, cCon2, initCon;
	// need to convert g/s/cm^3 to kg/s/m^3 by 1000.
	double conversion = 1000.;
	double x, xc, dx, x1, x2, s;
	int c1ID, c2ID, c3ID, c4ID, c5ID, c6ID;
	double c1Con, c2Con, c3Con, c4Con, c5Con, c6Con;
	double a1 = 8.30980E-04, a2 = 4.32710E-03, a3 = 4.19580E-03;
	double a4 = 1.19610E-02, a5 = 3.47340E-03, a6 = 1.22760E-03;
   double lambdaC1 = 0.0125, lambdaC2 = 0.0318, lambdaC3 = 0.109;
	double lambdaC4 = 0.3170, lambdaC5 = 1.3500, lambdaC6 = 8.640;
	std::vector<double> c1Coeffs = {-lambdaC1, 0.0, 0.0, 0.0, 0.0, 0.0};
	std::vector<double> c2Coeffs = {0.0, -lambdaC2, 0.0, 0.0, 0.0, 0.0};
	std::vector<double> c3Coeffs = {0.0, 0.0, -lambdaC3, 0.0, 0.0, 0.0};
	std::vector<double> c4Coeffs = {0.0, 0.0, 0.0, -lambdaC4, 0.0, 0.0};
	std::vector<double> c5Coeffs = {0.0, 0.0, 0.0, 0.0, -lambdaC5, 0.0};
	std::vector<double> c6Coeffs = {0.0, 0.0, 0.0, 0.0, 0.0, -lambdaC6};
	//std::vector<std::string> solvers {"CRAM", "parabolic", "hyperbolic"};
	//std::vector<std::string> solvers {"pade-method1", "pade-method2"};
	std::vector<std::string> solvers {"CRAM"};
	//std::vector<std::string> solvers {"BDF1", "BDF2", "BDF3", "BDF4", "BDF5", "BDF6"};
	meshCell* cell = nullptr;
	std::string outputFileName;
	
	// Build the Mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	// Sets the y velocity
	model.setConstantYVelocity(v);
	// Build species driver
	speciesDriver spec = speciesDriver(&model);
	// Loops over different solvers
	for (std::string &solverType : solvers){

		std::ofstream outputFile;
		outputFileName = "problem3"+solverType+".out";
		outputFile.open(outputFileName, std::ios::out | std::ios::trunc);

		// Loops over number of time steps
		for (double &numOfSteps	: steps){
			dt = tEnd/numOfSteps;
			outputFile << "TimeStepSize: "+std::to_string(dt)+"\n";

			// sets the solver
			spec.setMatrixExpSolver(solverType);
			//spec.setIntegratorSolver("implicit", solverType);
			// Add species
			c1ID = spec.addSpecies(1.0, 0.0, 0.0);
			c2ID = spec.addSpecies(1.0, 0.0, 0.0);
			c3ID = spec.addSpecies(1.0, 0.0, 0.0);
			c4ID = spec.addSpecies(1.0, 0.0, 0.0);
			c5ID = spec.addSpecies(1.0, 0.0, 0.0);
			c6ID = spec.addSpecies(1.0, 0.0, 0.0);

			// Set periodic BCs
			spec.setBoundaryCondition("periodic","north", c1ID);
			spec.setBoundaryCondition("periodic","south", c1ID);

			spec.setBoundaryCondition("periodic","north", c2ID);
			spec.setBoundaryCondition("periodic","south", c2ID);

			spec.setBoundaryCondition("periodic","north", c3ID);
			spec.setBoundaryCondition("periodic","south", c3ID);

			spec.setBoundaryCondition("periodic","north", c4ID);
			spec.setBoundaryCondition("periodic","south", c4ID);

			spec.setBoundaryCondition("periodic","north", c5ID);
			spec.setBoundaryCondition("periodic","south", c5ID);

			spec.setBoundaryCondition("periodic","north", c6ID);
			spec.setBoundaryCondition("periodic","south", c6ID);

			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);	

					// Calculates the x positions as the cell faces
					dx = cell->dy;
					xc = cell->y;
					x2 = xc + dx/2;
					x1 = xc - dx/2;

					// Calculates the initial concentration from MVT. 
					//initCon = (5./dx)*(x2*x2 - x1*x1);		
					initCon = 0.0;
					s = (1./dx)*(yLength/M_PI)*(cos(M_PI*x1/yLength) - 
						cos(M_PI*x2/yLength));

					spec.setSpeciesCon(i,j,c1ID, initCon);
					spec.setSpeciesCon(i,j,c2ID, initCon);
					spec.setSpeciesCon(i,j,c3ID, initCon);
					spec.setSpeciesCon(i,j,c4ID, initCon);
					spec.setSpeciesCon(i,j,c5ID, initCon);
					spec.setSpeciesCon(i,j,c6ID, initCon);

					// Sets the sourses
					// need to convert g/s/cm^3 to kg/s/m^3
					spec.setSpeciesSource(i, j, c1ID, c1Coeffs, s*a1*conversion);
					spec.setSpeciesSource(i, j, c2ID, c2Coeffs, s*a2*conversion);
					spec.setSpeciesSource(i, j, c3ID, c3Coeffs, s*a3*conversion);
					spec.setSpeciesSource(i, j, c4ID, c4Coeffs, s*a4*conversion);
					spec.setSpeciesSource(i, j, c5ID, c5Coeffs, s*a5*conversion);
					spec.setSpeciesSource(i, j, c6ID, c6Coeffs, s*a6*conversion);
				}
			}
			for (int step = 1; step <= numOfSteps; step++){
				t = step*dt;
				// Solve with CRAM
				//spec.solveImplicit(t);
				spec.solve(t);
			}
			// Get the solution
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);
					x = cell->y;
					c1Con = spec.getSpecies(i, j, c1ID);
					c2Con = spec.getSpecies(i, j, c2ID);
					c3Con = spec.getSpecies(i, j, c3ID);
					c4Con = spec.getSpecies(i, j, c4ID);
					c5Con = spec.getSpecies(i, j, c5ID);
					c6Con = spec.getSpecies(i, j, c6ID);
					outputFile << std::setprecision(16) << c1Con << " " << c2Con << " " << c3Con 
					//std::cout << i << " " << j << " " << c1Con << " " << c2Con << " " << c3Con
					<< " " << c4Con << " " << c5Con << " " << c6Con << std::endl;
				}
			}
		spec.clean();
		}
	spec.clean();
	}
}
//*****************************************************************************
// Test that the species driver sets up the problem right and solves the 
// system right. 
//
//	Note: The units for this problem kinda dont matter either. This is just an
//	internal ODE test, the ODE solution is provided and calculated using the 
//	same problem units. 
//*****************************************************************************
void testXenonIodineNoFlow(int myid){
	int xCells = 1, yCells = 1;
	double xLength = 1.0, yLength = 1.0;
	double xenonInitCon = 0.0, iodineInitCon = 0.0;
	double xenonMM = 135.0, iodineMM = 135.0;
	double numOfSteps = 4.;
	double tEnd = 1000.0;
	double t;
	double dt = tEnd/numOfSteps;
   double lambda_I = 2.11E-5;
   double lambda_xe = 2.9306E-5;
   double sigma_a = 2.002E-22;
   double Sigma_f = 9.7532E-1;
   double flux = 2.5E16;
   double gamma_xe = 0.002468;
   double gamma_I = 0.063033;
	double D_xe = 0.0, D_I = 0.0;
	double N_xe_0 = 0.0, N_I_0 = 0.0;
	int xenonID, iodineID;
	double xenonCon, iodineCon;
	std::vector<double> xenonCoeffs = {-lambda_xe-sigma_a*flux, lambda_I};
	std::vector<double> iodineCoeffs = {0.0, -lambda_I};
	double xenonS = gamma_xe*Sigma_f*flux;
	double iodineS = gamma_I*Sigma_f*flux;

	modelMesh model(xCells, yCells, xLength, yLength);
	speciesDriver spec = speciesDriver(&model);

	xenonID = spec.addSpecies(xenonMM, N_xe_0, D_xe);
	iodineID = spec.addSpecies(iodineMM, N_I_0, D_I);
	//spec.setIntegratorSolver("explicit", "classic fourth-order");

	// Set source
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			spec.setSpeciesSource(i, j, xenonID, xenonCoeffs, xenonS);
			spec.setSpeciesSource(i, j, iodineID, iodineCoeffs, iodineS);
		}
	}
	
	for (int step = 1; step <= numOfSteps; step++){
		t = step*dt;
		// Solve with CRAM
		spec.solve(t);

		// Stuff for analytical solution
		double a = lambda_xe + sigma_a*flux;
   	double b = gamma_I*Sigma_f*flux;
   	double d = lambda_I*N_I_0;
   	double k = N_xe_0 - (d-b)/(a - lambda_I) - (b + gamma_xe*Sigma_f*flux)/a;

   	// Xenon solution
   	double N_xe = -b/(a-lambda_I)*exp(-lambda_I*t) + b/a +
   	   d*exp(-lambda_I*t)/(a - lambda_I) + k*exp(-a*t) +
   	   gamma_xe*Sigma_f*flux/a;

   	// Iodine solution
   	double N_I = b/lambda_I*(1. - exp(-lambda_I*t)) + N_I_0*exp(-lambda_I*t);

		// Gets species Concentrations
		if (myid==0){
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					xenonCon = spec.getSpecies(i, j, xenonID);
					iodineCon = spec.getSpecies(i, j, iodineID);
					assert(isApprox(xenonCon, N_xe, 1.e5, 1.e-11));
					assert(isApprox(iodineCon, N_I, 1.e5, 1.e-11));
					//std::cout << xenonCon << std::endl; 
					//std::cout << iodineCon << std::endl; 
				}
			}
		}
	}
	
	model.clean();
	spec.clean();
}
//*****************************************************************************
// Test Xenon iodine flow problem in the y direction
//
//	Note: The units for this problem kinda dont matter either. This is just an
//	internal ODE test, the ODE solution is provided and calculated using the 
//	same problem units. 
//*****************************************************************************
void testXenonIodineYFlow(int myid){
	int xCells = 1, yCells = 500;
	double xLength = 0.0, yLength = 10.0;
	double yVelocity = 8.0;
	double xenonInitCon = 5e-6, iodineInitCon = 5e-6;
	double t = 10000000.0;
	double xenonMM = 135.0, iodineMM = 135.0;
	double AvogNum = 6.02214076E23;
   double lambda_I = 2.11E-5;
   double lambda_xe = 2.9306E-5;
   double sigma_a = 2.002E-22;
   double Sigma_f = 9.7532E-1;
   double flux = 2.5E16;
   double gamma_xe = 0.002468;
   double gamma_I = 0.063033;
	double D_xe = 0.0, D_I = 0.0;
	double N_xe_0 = 0.0, N_I_0 = 0.0;
	int xenonID, iodineID;
	double xenonCon, iodineCon, error = 0.0;
	std::vector<double> xenonCoeffs = {-lambda_xe-sigma_a*flux, lambda_I};
	std::vector<double> iodineCoeffs = {0.0, -lambda_I};
	double xenonS = gamma_xe*Sigma_f*flux*xenonMM/AvogNum;
	double iodineS = gamma_I*Sigma_f*flux*iodineMM/AvogNum;

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	model.addBoundarySurface("south");
	model.addBoundarySurface("north");

	// Sets the x velocity
	model.setConstantYVelocity(yVelocity);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds xenon and iodine species
	xenonID = spec.addSpecies(xenonMM, N_xe_0, D_xe);
	iodineID = spec.addSpecies(iodineMM, N_I_0, D_I);
	spec.setBoundaryCondition("dirichlet", "south", xenonID, xenonInitCon);
	spec.setBoundaryCondition("dirichlet","south", iodineID, iodineInitCon);
	spec.setBoundaryCondition("newmann", "north", xenonID, 0.0);
	spec.setBoundaryCondition("newmann","north", iodineID, 0.0);

	// Set source
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			spec.setSpeciesSource(i, j, xenonID, xenonCoeffs, xenonS);
			spec.setSpeciesSource(i, j, iodineID, iodineCoeffs, iodineS);
		}
	}

	spec.solve(t);

	// Gets species Concentrations
	if (myid==0){
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				xenonCon = spec.getSpecies(i, j, xenonID);
				iodineCon = spec.getSpecies(i, j, iodineID);
				meshCell* cell = model.getCellByLoc(i,j);
				double y = cell->y;
				// Iodine solution
				double b = gamma_I*Sigma_f*flux*iodineMM/AvogNum/lambda_I;
   			double N_I = b + (iodineInitCon - b)*exp(-lambda_I/yVelocity*y);

				//std::cout << y << " " << iodineCon << " " << N_I << std::endl;
				error = std::max(std::abs(iodineCon - N_I)/N_I, error);
				//error = std::abs(iodineCon - N_I)/N_I;
				//std::cout << y << " " << error << std::endl;
				assert(isApprox(iodineCon, N_I, 1.e-9, 1.e-4));
			}
		}
	}
	//std::cout << "Max l-1 error: " << error << std::endl;

	
	model.clean();
	spec.clean();

}

//*****************************************************************************
// Test Xenon iodine flow problem in the x direction
//
//	Note: The units for this problem kinda dont matter either. This is just an
//	internal ODE test, the ODE solution is provided and calculated using the 
//	same problem units. 
//*****************************************************************************
void testXenonIodineXFlow(int myid){
	int xCells = 500, yCells = 1;
	double xLength = 10.0, yLength = 0.0;
	double xVelocity = 8.0;
	double xenonInitCon = 5e-6, iodineInitCon = 5e-6;
	double xenonMM = 135.0, iodineMM = 135.0;
	double AvogNum = 6.02214076E23;
	double t = 10000000.0;
   double lambda_I = 2.11E-5;
   double lambda_xe = 2.9306E-5;
   double sigma_a = 2.002E-22;
   double Sigma_f = 9.7532E-1;
   double flux = 2.5E16;
   double gamma_xe = 0.002468;
   double gamma_I = 0.063033;
	double D_xe = 0.0, D_I = 0.0;
	double N_xe_0 = 0.0, N_I_0 = 0.0;
	int xenonID, iodineID;
	double xenonCon, iodineCon, error = 0.0;
	double xenonError, iodineError;
	std::vector<double> xenonCoeffs = {-lambda_xe-sigma_a*flux, lambda_I};
	std::vector<double> iodineCoeffs = {0.0, -lambda_I};
	double xenonS = gamma_xe*Sigma_f*flux*xenonMM/AvogNum;
	double iodineS = gamma_I*Sigma_f*flux*iodineMM/AvogNum;

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	model.addBoundarySurface("west");
	model.addBoundarySurface("east");

	// Sets the x velocity
	model.setConstantXVelocity(xVelocity);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds xenon and iodine species
	xenonID = spec.addSpecies(xenonMM, N_xe_0, D_xe);
	iodineID = spec.addSpecies(iodineMM, N_I_0, D_I);
	spec.setBoundaryCondition("dirichlet","west", xenonID, xenonInitCon);
	spec.setBoundaryCondition("dirichlet","west", iodineID, iodineInitCon);
	spec.setBoundaryCondition("newmann","east", xenonID, 0.0);
	spec.setBoundaryCondition("newmann","east", iodineID, 0.0);

	// Set source
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			spec.setSpeciesSource(i, j, xenonID, xenonCoeffs, xenonS);
			spec.setSpeciesSource(i, j, iodineID, iodineCoeffs, iodineS);
		}
	}

	// Solve with CRAM
	spec.solve(t);

	// Gets species Concentrations
	if (myid==0){
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				xenonCon = spec.getSpecies(i, j, xenonID);
				iodineCon = spec.getSpecies(i, j, iodineID);
				meshCell* cell = model.getCellByLoc(i,j);
				double x = cell->x;
				// Iodine solution
				double b = gamma_I*Sigma_f*flux*iodineMM/AvogNum/lambda_I;
   			double N_I = b + (iodineInitCon - b)*exp(-lambda_I/xVelocity*x);

				//std::cout << xenonCon << " " << iodineCon << " " << N_I << std::endl;
				//iodineError = std::abs(iodineCon-N_I)/N_I;
				//printf (" %2i %2i %e \n", i, j, iodineError);
				//error = std::max(std::abs(iodineCon - N_I)/N_I, error);
				//std::cout << x << " " << iodineError << std::endl;
				//std::cout << x << " " << iodineCon << " " << N_I << std::endl;
				assert(isApprox(iodineCon, N_I, 1.e-9, 1.e-4));
			}
		}
	}
	//std::cout << "Max l-1 error: " << error << std::endl;

	
	model.clean();
	spec.clean();

}
//*****************************************************************************
// Test 2D diffusion
//
//	Note: The units for this problem kinda dont matter either. This is just an
//	internal ODE test, the ODE solution is provided and calculated using the 
//	same problem units. 
//*****************************************************************************
void testDiffusion1(int myid){
	int xCells = 50, yCells = 50;
	double xLength = 15.0, yLength = 20.0;
	double totalTime = 1000000.0;
	double t = 0.0;
	int steps = 1;
	double dt = totalTime/steps;
	double D_spec = 1.0;
	int specID;
	double specCon;
	double s, y, x;
	double exact, temp1, temp2, temp3, error;

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	model.addBoundarySurface("north");
	model.addBoundarySurface("south");
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds xenon and iodine species
	specID = spec.addSpecies(1.0, 0.0, D_spec);
	spec.setBoundaryCondition("dirichlet","south", specID, 0.0);
	spec.setBoundaryCondition("dirichlet","east", specID, 0.0);
	spec.setBoundaryCondition("dirichlet","west", specID, 0.0);
	spec.setBoundaryCondition("dirichlet","north", specID, 100.0);

	// Set source
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			//s = exp(-20.*pow(x-1.0/2.,2) - 20.*pow(y-1.0/2.,2));
			s = 0.0;
			spec.setSpeciesCon(i, j, specID, s);
		}
	}

	error = 0.0;
	t = 0.0;
	for (int k = 0; k < steps; k++){
		t = t + dt;
		spec.solve(t);

		std::ofstream outputFile;
		outputFile.open("Diffusion1.out", std::ios::out | std::ios::trunc);
		outputFile << "Time: "+std::to_string(t)+"\n";

		// Gets species Concentrations
		if (myid==0){
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					specCon = spec.getSpecies(i, j, specID);
					meshCell* cell = model.getCellByLoc(i,j);
					x = cell->x;
					y = cell->y;
					exact = 0.0;

					// Calculates exact solution from the first project
					// of my CFD class
					for (int k = 1; k < 102; k += 2){
						temp1 = sinh(k*M_PI*y/xLength);
						temp2 = sin(k*M_PI*x/xLength);
						temp3 = 1./(k*sinh(k*M_PI*yLength/xLength));
						exact += temp1*temp2*temp3;						
					}
					exact = exact*400./M_PI;
					error += pow(std::abs(specCon - exact), 2);
					outputFile << i << " " << j << " " << specCon << std::endl;
				}
			}
		}
		error = pow(error,0.5)/(xCells*yCells);
		assert(error < 0.003);
	}

	
	model.clean();
	spec.clean();
	
}
//*****************************************************************************
// Test 2D diffusion
//
//*****************************************************************************
void testDiffusion2(int myid){
	int xCells = 50, yCells = 50;
	double xLength = 1.0, yLength = 1.0;
	double totalTime = 2.0;
	double t = 0.0;
	int steps = 1;
	double dt = totalTime/steps;
	double D_spec = 1./(4.*2.*M_PI*M_PI);
	int specID;
	double specCon;
	double s, sx, sy, yc, y1, y2, xc, x1, x2, dx, dy;
	double exact, temp1, temp2, temp3, error;

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Sets the species matrix exp solver
	//spec.setMatrixExpSolver("taylor");

	// Adds xenon and iodine species
	specID = spec.addSpecies(1.0, 0.0, D_spec);
	spec.setBoundaryCondition("periodic","south", specID);
	spec.setBoundaryCondition("periodic","east", specID);
	spec.setBoundaryCondition("periodic","west", specID);
	spec.setBoundaryCondition("periodic","north", specID);

	// Set inital condition
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			meshCell* cell = model.getCellByLoc(i,j);
			dx = cell->dx, dy = cell->dy;
			xc = cell->x, yc = cell->y;
			x2 = xc + dx/2, x1 = xc - dx/2;
			y2 = yc + dy/2, y1 = yc - dy/2;
			sx = 1./(2.*M_PI)*(cos(2.*M_PI*x1) - cos(2.*M_PI*x2));
			sy = 1./(2.*M_PI)*(cos(2.*M_PI*y1) - cos(2.*M_PI*y2));
			s = (1./dx)*(1./dy)*sx*sy;
			spec.setSpeciesCon(i, j, specID, s);
		}
	}

	error = 0.0;
	t = 0.0;
	for (int k = 0; k < steps; k++){
		t = t + dt;
		// Solve
		spec.solve(t);

		std::ofstream outputFile;
		outputFile.open("Diffusion2.out", std::ios::out | std::ios::trunc);
		outputFile << "Time: "+std::to_string(t)+"\n";

		// Gets species Concentrations
		if (myid==0){
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					specCon = spec.getSpecies(i, j, specID);
					meshCell* cell = model.getCellByLoc(i,j);
					xc = cell->x;
					yc = cell->y;
					exact = exp(-t)*sin(2.*M_PI*xc)*sin(2.*M_PI*yc);

					error += pow(std::abs(specCon - exact), 2);
					//std::cout << i << " " << j << " " << specCon << " " << exact << std::endl;
				}
			}
		}
		error = error/(xCells*yCells);
		assert(error < 1.0e-8);
	}

	
	model.clean();
	spec.clean();
	
}
//*****************************************************************************
// Test gas sparging model with no flow
//
// dC_gas/dt = [K*a/V* (C*_a - C_liq/(1-alpha)) ]
// dC_liq/dt = [K*a/V* (C_liq/(1-alpha) - C*_a) ]
//*****************************************************************************
void testGasTransportNoFlow(int myid){
	int xCells = 10, yCells = 1;
	double xLength = 500.0, yLength = 0.0;
	double totalTime = 100.0; // seconds
	double t = 0.0;
	int steps = 1;
	double dt = totalTime/steps;
	std::string limiter = "First order upwind", solverType = "hyperbolic";
	int liq1ID, gas1ID, liq2ID, gas2ID;
	double liq1Con, gas1Con, liq2Con, gas2Con, x;
	double temperature = 600.0;		// kelvin
	double gasIntAreaCon = 50.0;		// 1/m
	double gasVoidFraction = 0.001;	// fraction
	double pressure = 101325.0;		// Pa
	double k = 3.0;						// m/s
	double H = 2.71e-8;					// mol/m^3/Pa
	double spec1InitCon = 10., spec2InitCon = 5.;

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");

	// Sets scalar data
	model.setSystemTemperature(temperature);
	model.setSystemPressure(pressure);
	model.setSystemGasInterfacialAreaCon(gasIntAreaCon);
	model.setSystemGasVoidFraction(gasVoidFraction);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Sets the flux limiter type
	spec.setFluxLimiter(limiter);
	// Sets the Matrix expential solver
	spec.setMatrixExpSolver(solverType);

	// Adds gas and liquid spec
	liq1ID = spec.addSpecies(1.0, spec1InitCon, 0.0);
	liq2ID = spec.addSpecies(1.0, spec2InitCon, 0.0);
	gas1ID = spec.addSpecies(1.0, 0.0, 0.0);
	gas2ID = spec.addSpecies(1.0, 0.0, 0.0);

	// Adds gas sparging model
	spec.setGasSparging({k, 0.5*k}, {H, 1.e-2*H}, {liq1ID, liq2ID}, {gas1ID, gas2ID});

	for (int k = 0; k < steps; k++){
		t = t + dt;
		// Solve with CRAM
		spec.solve(t);
	}

	// Gets species Concentrations
	if (myid==0){
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				liq1Con = spec.getSpecies(i, j, liq1ID);
				gas1Con = spec.getSpecies(i, j, gas1ID);
				liq2Con = spec.getSpecies(i, j, liq2ID);
				gas2Con = spec.getSpecies(i, j, gas2ID);
				meshCell* cell = model.getCellByLoc(i,j);
				x = cell->x;
				assert(isApprox(spec1InitCon, liq1Con + gas1Con));
				assert(isApprox(spec2InitCon, liq2Con + gas2Con));
			}
		}
	}
	
	model.clean();
	spec.clean();

}

//*****************************************************************************
// Test gas sparging model with flow
//
// dC_gas/dt = (1/v_x)*[K*a/V* (C*_a - C_liq/(1-alpha)) ]
// dC_liq/dt = (1/v_x)*[K*a/V* (C_liq/(1-alpha) - C*_a) ]
//*****************************************************************************
void testGasTransportFlow(int myid){
	int xCells = 10, yCells = 1;
	double xLength = 500.0, yLength = 0.0;
	double totalTime = 100.0; // seconds
	double t = 0.0;
	int steps = 1;
	double dt = totalTime/steps;
	std::string limiter = "First order upwind", solverType = "hyperbolic";
	int liq1ID, gas1ID, liq2ID, gas2ID;
	double liq1Con, gas1Con, liq2Con, gas2Con, x;
	double temperature = 600.0;		// kelvin
	double gasIntAreaCon = 50.0;		// 1/m
	double gasVoidFraction = 0.001;	// fraction
	double pressure = 101325.0;		// Pa
	double k = 3.0;						// m/s
	double H = 2.71e-8;					// mol/m^3/Pa
	double v_x = 1.5;						// m/s

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");

	// Sets scalar data
	model.setSystemTemperature(temperature);
	model.setSystemPressure(pressure);
	model.setSystemGasInterfacialAreaCon(gasIntAreaCon);
	model.setSystemGasVoidFraction(gasVoidFraction);
	model.setConstantXVelocity(v_x);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Sets the flux limiter type
	spec.setFluxLimiter(limiter);
	// Sets the Matrix expential solver
	spec.setMatrixExpSolver(solverType);

	// Adds gas and liquid spec
	liq1ID = spec.addSpecies(1.0, 0.0, 0.0);
	liq2ID = spec.addSpecies(1.0, 0.0, 0.0);
	gas1ID = spec.addSpecies(1.0, 0.0, 0.0);
	gas2ID = spec.addSpecies(1.0, 0.0, 0.0);

	// Dirichlet Bc as the pipe inlet. Outflow Bc as the pipe exit
	spec.setBoundaryCondition("dirichlet","west", {liq1ID, liq2ID, gas1ID, gas2ID}, {10.0, 5.0, 0.0, 0.0});
	spec.setBoundaryCondition("newmann","east", {liq1ID, liq2ID, gas1ID, gas2ID}, {0.0, 0.0, 0.0, 0.0});

	// Adds gas sparging model
	spec.setGasSparging({k, 0.5*k}, {H, 1.e-2*H}, {liq1ID, liq2ID}, {gas1ID, gas2ID});

	for (int k = 0; k < steps; k++){
		t = t + dt;
		// Solve with CRAM
		spec.solve(t);
	}

	// Gets species Concentrations
	if (myid==0){
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				liq1Con = spec.getSpecies(i, j, liq1ID);
				gas1Con = spec.getSpecies(i, j, gas1ID);
				liq2Con = spec.getSpecies(i, j, liq2ID);
				gas2Con = spec.getSpecies(i, j, gas2ID);
				meshCell* cell = model.getCellByLoc(i,j);
				x = cell->x;
				printf("%i %i %f %f %f %f\n", i, j, liq1Con, gas1Con, liq2Con, gas2Con);
				//assert(liq1Con < gas1Con);
				//assert(liq2Con < gas2Con);
			}
		}
	}
	
	model.clean();
	spec.clean();

}

//*****************************************************************************
// Main test
//*****************************************************************************
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	//testProblem1(myid);
	//testProblem2(myid);
	//testProblem2Krylov(myid);
	//testProblem2IntegratorMethods(myid);
	//testProblem3(myid);
	//testXenonIodineNoFlow(myid);
	//testXenonIodineYFlow(myid);
	//testXenonIodineXFlow(myid);
	//testDiffusion1(myid);
	//testDiffusion2(myid);
	testGasTransportNoFlow(myid);
	//testGasTransportFlow(myid);

	mpi.finalize();
}
