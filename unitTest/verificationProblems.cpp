#define _USE_MATH_DEFINES
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>
#include <assert.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>

#include "mpiProcess.h"
#include "modelMesh.h"
#include "speciesDriver.h"
#include "meshCellData.h"
#include "species.h"

using namespace Eigen;
//*****************************************************************************
// Test if two number are approx equal
//
// @param goalVal		The real solution value
// @param testVal		Test value
// @param rtol			Relative tolerance
// @param atol			Absolution tolerance
//
// Values for rtol and atol were taken from the default values for numpys 
// isApprox function.
//*****************************************************************************
bool isApprox(double goalVal, double testVal, double rtol = 1e-5, double atol = 1e-8){
	bool retBool = false;
	bool rtolBool = false;
	bool atolBool = false;

	double diff = abs(goalVal - testVal);
	if (diff < rtol) {rtolBool = true;};
	if (diff/goalVal < atol) {atolBool = true;};
	if (rtolBool and atolBool) {retBool = true;};
	return retBool;
}

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
//*****************************************************************************
void testProblem1NoFlow(int myid){
	int xCells = 1, yCells = 1;
	double xLength = 1.0, yLength = 1.0;
	double N1InitCon = 10000.0, N2InitCon = 0.0, N3InitCon = 0.0;
	double N1MM = 1.0, N2MM = 1.0, N3MM = 1.0;
	double numOfSteps = 100;
	double tEnd = 850.0;
	double t;
	double dt = tEnd/numOfSteps;
   double lambda1 = 1.0/1.0e2, lambda2 = 0.5/1.0e2, lambda3 = 1.1/1.0e2;
	double DN1 = 0.0, DN2 = 0.0, DN3 = 0.0;
	int N1ID, N2ID, N3ID;
	double N1Con, N2Con, N3Con;
	std::vector<double> N1Coeffs = {-lambda1, 0.0, lambda3};
	std::vector<double> N2Coeffs = {lambda1, -lambda2, 0.0};
	std::vector<double> N3Coeffs = {0.0, lambda2, -lambda3};

	modelMesh model(xCells, yCells, xLength, yLength);
	speciesDriver spec = speciesDriver(&model);
	N1ID = spec.addSpecies(N1MM, N1InitCon, DN1);
	N2ID = spec.addSpecies(N2MM, N2InitCon, DN2);
	N3ID = spec.addSpecies(N3MM, N3InitCon, DN3);

	// Set source
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			spec.setSpeciesSource(i, j, N1ID, N1Coeffs, 0.0);
			spec.setSpeciesSource(i, j, N2ID, N2Coeffs, 0.0);
			spec.setSpeciesSource(i, j, N3ID, N3Coeffs, 0.0);
		}
	}
	
	for (int step = 1; step <= numOfSteps; step++){
		t = step*dt;
		// Solve with CRAM
		spec.solve(t);

		std::ofstream outputFile;
		outputFile.open("problem1.out", std::ios_base::app);
		//outputFile << "Time: "+std::to_string(t)+"\n";
		//std::cout.precision(16);

		// Gets species Concentrations
		if (myid==0){
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					N1Con = spec.getSpecies(i, j, N1ID);
					N2Con = spec.getSpecies(i, j, N2ID);
					N3Con = spec.getSpecies(i, j, N3ID);
					//std::cout << t << " " << N1Con << std::endl;
					//outputFile <<  std::setprecision(16) << t << " " << N1Con 
					//<< " " << N2Con << " " << N3Con <<  std::endl;
				}
			}
		}
	}
	
	model.clean();
	spec.clean();
}
//*****************************************************************************
// Test that the species driver sets up the problem right and solves the 
// system right. 
//*****************************************************************************
void testXenonIodineNoFlow(int myid){
	int xCells = 1, yCells = 1;
	double xLength = 1.0, yLength = 1.0;
	double xenonInitCon = 0.0, iodineInitCon = 0.0;
	double xenonMM = 135.0, iodineMM = 135.0;
	double numOfSteps = 10.;
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
					std::cout << std::abs(iodineCon - N_I)/N_I << std::endl;
					std::cout << std::abs(xenonCon - N_xe)/N_xe <<  std::endl;
					assert(isApprox(xenonCon, N_xe, 1.e5, 1.e-11));
					assert(isApprox(iodineCon, N_I, 1.e5, 1.e-11));
				}
			}
		}
	}
	
	model.clean();
	spec.clean();
}
//*****************************************************************************
// Test Xenon iodine flow problem in the y direction
//*****************************************************************************
void testXenonIodineYFlow(int myid){
	int xCells = 1, yCells = 500;
	double xLength = 1.0, yLength = 10.0;
	double yVelocity = 8.0;
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
	std::vector<double> xenonCoeffs = {-lambda_xe-sigma_a*flux, lambda_I};
	std::vector<double> iodineCoeffs = {0.0, -lambda_I};
	double xenonS = gamma_xe*Sigma_f*flux*xenonMM/AvogNum;
	double iodineS = gamma_I*Sigma_f*flux*iodineMM/AvogNum;

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);

	// Sets the x velocity
	model.setConstantYVelocity(yVelocity);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds xenon and iodine species
	xenonID = spec.addSpecies(xenonMM, N_xe_0, D_xe);
	iodineID = spec.addSpecies(iodineMM, N_I_0, D_I);
	spec.setBoundaryCondition("dirichlet", "south", xenonID, xenonInitCon);
	spec.setBoundaryCondition("dirichlet","south", iodineID, iodineInitCon);

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
				double y = cell->y;
				// Iodine solution
				double b = gamma_I*Sigma_f*flux*iodineMM/AvogNum/lambda_I;
   			double N_I = b + (iodineInitCon - b)*exp(-lambda_I/yVelocity*y);

				//std::cout << xenonCon << " " << iodineCon << " " << N_I << std::endl;
				//error = std::max(std::abs(iodineCon - N_I)/N_I, error);
				//std::cout << y << " " << error << std::endl;
				assert(isApprox(iodineCon, N_I));
			}
		}
	}
	//std::cout << "Max l-1 error: " << error << std::endl;

	
	model.clean();
	spec.clean();

}

//*****************************************************************************
// Test Xenon iodine flow problem in the x direction
//*****************************************************************************
void testXenonIodineXFlow(int myid){
	int xCells = 500, yCells = 1;
	double xLength = 10.0, yLength = 1.0;
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

	// Sets the x velocity
	model.setConstantXVelocity(xVelocity);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds xenon and iodine species
	xenonID = spec.addSpecies(xenonMM, N_xe_0, D_xe);
	iodineID = spec.addSpecies(iodineMM, N_I_0, D_I);
	spec.setBoundaryCondition("dirichlet","west", xenonID, xenonInitCon);
	spec.setBoundaryCondition("dirichlet","west", iodineID, iodineInitCon);

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
				assert(isApprox(iodineCon, N_I));
			}
		}
	}
	//std::cout << "Max l-1 error: " << error << std::endl;

	
	model.clean();
	spec.clean();

}
//*****************************************************************************
// Test 2D diffusion
//*****************************************************************************
void testDiffusion2D(int myid){
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
	for (int k = 0; k < steps; k++){
		t = t + dt;
		//std::cout << t << std::endl;
		// Solve with CRAM
		spec.solve(t);

		std::ofstream outputFile;
		outputFile.open("Diffusion2D.out", std::ios_base::app);
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
					//std::cout << x << " " << y << " " << specCon << " " << exact << std::endl;
					outputFile << i << " " << j << " " << specCon << std::endl;
				}
			}
		}
		error = pow(error,0.5)/(xCells*yCells);
		assert(error < 0.008);
		//std::cout << "Error: " << error << std::endl;
	}

	
	model.clean();
	spec.clean();

}
//*****************************************************************************
// Test neutron precursors for sinlge channel
//*****************************************************************************
void testNeutronPrecursorsFlow(int myid){
	double t = 0.0;
	int steps = 1;
	double totalTime = 50.0;
	double dt = totalTime/steps;
	int xCells = 1, yCells = 100;
	double xLength = 1.0, yLength = 15.5;
	double scale;
	double AvogNum = 6.02214076E23;
	int c1ID, c2ID, c3ID, c4ID, c5ID, c6ID;
	double c1Con, c2Con, c3Con, c4Con, c5Con, c6Con;
	double c1Ana, c2Ana, c3Ana, c4Ana, c5Ana, c6Ana;
   double lambdaC1 = 0.0125, lambdaC2 = 0.0318, lambdaC3 = 0.109;
	double lambdaC4 = 0.3170, lambdaC5 = 1.3500, lambdaC6 = 8.640;
	double c1InitCon = 0.0, c2InitCon = 0.0, c3InitCon = 0.0, c4InitCon = 0.0;
	double c5InitCon = 0.0, c6InitCon = 0.0;
	double D_c1 = 0.0, D_c2 = 0.0, D_c3 = 0.0, D_c4 = 0.0, D_c5 = 0.0, D_c6 = 0.0;
	double y, s, y1, y2;
	double a1 = 8.30980E-04, a2 = 4.32710E-03, a3 = 4.19580E-03;
	double a4 = 1.19610E-02, a5 = 3.47340E-03, a6 = 1.22760E-03;
	double c1Error, c2Error, c3Error, c4Error, c5Error, c6Error;
	double c1Errorl1=0.0, c2Errorl1=0.0, c3Errorl1=0.0, c4Errorl1=0.0, c5Errorl1=0.0, c6Errorl1=0.0;
	double yVelocity = 7.0;
   MatrixD coeff(16,7);
	std::ofstream outputFile;
	std::vector<double> c1Coeffs = {-lambdaC1, 0.0, 0.0, 0.0, 0.0, 0.0};
	std::vector<double> c2Coeffs = {0.0, -lambdaC2, 0.0, 0.0, 0.0, 0.0};
	std::vector<double> c3Coeffs = {0.0, 0.0, -lambdaC3, 0.0, 0.0, 0.0};
	std::vector<double> c4Coeffs = {0.0, 0.0, 0.0, -lambdaC4, 0.0, 0.0};
	std::vector<double> c5Coeffs = {0.0, 0.0, 0.0, 0.0, -lambdaC5, 0.0};
	std::vector<double> c6Coeffs = {0.0, 0.0, 0.0, 0.0, 0.0, -lambdaC6};
	//std::vector<double> c5Coeffs = { -lambdaC5, 0.0};
	//std::vector<double> c6Coeffs = { 0.0, -lambdaC6};

   // Coefficients in order of precursor groups for columns
	coeff(0,0)  = 1.9490E-04, coeff(0,1)  = 1.0149E-03, coeff(0,2)  = 9.8409E-04;
	coeff(1,0)  = 3.5140E-04, coeff(1,1)  = 1.8298E-03, coeff(1,2)  = 1.7743E-03;
	coeff(2,0)  = 4.8659E-04, coeff(2,1)  = 2.5338E-03, coeff(2,2)  = 2.4569E-03;
	coeff(3,0)  = 6.0108E-04, coeff(3,1)  = 3.1300E-03, coeff(3,2)  = 3.0350E-03;
	coeff(4,0)  = 6.9421E-04, coeff(4,1)  = 3.6150E-03, coeff(4,2)  = 3.5052E-03;
	coeff(5,0)  = 7.6455E-04, coeff(5,1)  = 3.9812E-03, coeff(5,2)  = 3.8604E-03;
	coeff(6,0)  = 8.1053E-04, coeff(6,1)  = 4.2207E-03, coeff(6,2)  = 4.0926E-03;
	coeff(7,0)  = 8.3098E-04, coeff(7,1)  = 4.3271E-03, coeff(7,2)  = 4.1958E-03;
	coeff(8,0)  = 8.2533E-04, coeff(8,1)  = 4.2977E-03, coeff(8,2)  = 4.1673E-03;
	coeff(9,0)  = 7.9376E-04, coeff(9,1)  = 4.1333E-03, coeff(9,2)  = 4.0079E-03;
	coeff(10,0) = 7.3720E-04, coeff(10,1) = 3.8388E-03, coeff(10,2) = 3.7223E-03;
	coeff(11,0) = 6.5731E-04, coeff(11,1) = 3.4228E-03, coeff(11,2) = 3.3189E-03;
	coeff(12,0) = 5.5629E-04, coeff(12,1) = 2.8967E-03, coeff(12,2) = 2.8088E-03;
	coeff(13,0) = 4.3661E-04, coeff(13,1) = 2.2735E-03, coeff(13,2) = 2.2045E-03;
	coeff(14,0) = 3.0078E-04, coeff(14,1) = 1.5662E-03, coeff(14,2) = 1.5187E-03;
	coeff(15,0) = 1.5123E-04, coeff(15,1) = 7.8752E-04, coeff(15,2) = 7.6362E-04;

   coeff(0,3)  = 2.8054E-03, coeff(0,4)  = 8.1466E-04, coeff(0,5)  = 2.8792E-04;
	coeff(1,3)  = 5.0582E-03, coeff(1,4)  = 1.4688E-03, coeff(1,5)  = 5.1911E-04;
	coeff(2,3)  = 7.0042E-03, coeff(2,4)  = 2.0339E-03, coeff(2,5)  = 7.1883E-04;
	coeff(3,3)  = 8.6521E-03, coeff(3,4)  = 2.5124E-03, coeff(3,5)  = 8.8795E-04;
	coeff(4,3)  = 9.9927E-03, coeff(4,4)  = 2.9017E-03, coeff(4,5)  = 1.0255E-03;
	coeff(5,3)  = 1.1005E-02, coeff(5,4)  = 3.1957E-03, coeff(5,5)  = 1.1294E-03;
	coeff(6,3)  = 1.1667E-02, coeff(6,4)  = 3.3880E-03, coeff(6,5)  = 1.1974E-03;
	coeff(7,3)  = 1.1961E-02, coeff(7,4)  = 3.4734E-03, coeff(7,5)  = 1.2276E-03;
	coeff(8,3)  = 1.1880E-02, coeff(8,4)  = 3.4498E-03, coeff(8,5)  = 1.2192E-03;
	coeff(9,3)  = 1.1426E-02, coeff(9,4)  = 3.3178E-03, coeff(9,5)  = 1.1726E-03;
	coeff(10,3) = 1.0611E-02, coeff(10,4) = 3.0814E-03, coeff(10,5) = 1.0890E-03;
	coeff(11,3) = 9.4616E-03, coeff(11,4) = 2.7475E-03, coeff(11,5) = 9.7103E-04;
	coeff(12,3) = 8.0074E-03, coeff(12,4) = 2.3252E-03, coeff(12,5) = 8.2178E-04;
	coeff(13,3) = 6.2846E-03, coeff(13,4) = 1.8250E-03, coeff(13,5) = 6.4498E-04;
	coeff(14,3) = 4.3295E-03, coeff(14,4) = 1.2572E-03, coeff(14,5) = 4.4433E-04;
	coeff(15,3) = 2.1769E-03, coeff(15,4) = 6.3215E-04, coeff(15,5) = 2.2341E-04;

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);

	// Sets the y velocity
	model.setConstantYVelocity(yVelocity);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds species
	c1ID = spec.addSpecies(1.0, 0.0, D_c1);
	c2ID = spec.addSpecies(1.0, 0.0, D_c2);
	c3ID = spec.addSpecies(1.0, 0.0, D_c3);
	c4ID = spec.addSpecies(1.0, 0.0, D_c4);
	c5ID = spec.addSpecies(1.0, 0.0, D_c5);
	c6ID = spec.addSpecies(1.0, 0.0, D_c6);

	// Sets BCs
	spec.setBoundaryCondition("dirichlet", "south", c1ID, c1InitCon);
	spec.setBoundaryCondition("dirichlet", "south", c2ID, c2InitCon);
	spec.setBoundaryCondition("dirichlet", "south", c3ID, c3InitCon);
	spec.setBoundaryCondition("dirichlet", "south", c4ID, c4InitCon);
	spec.setBoundaryCondition("dirichlet", "south", c5ID, c5InitCon);
	spec.setBoundaryCondition("dirichlet", "south", c6ID, c6InitCon);
	// Sets BCs
	//spec.setBoundaryCondition("periodic","south", c1ID, c1InitCon);
	//spec.setBoundaryCondition("periodic","south", c2ID, c2InitCon);
	//spec.setBoundaryCondition("periodic","south", c3ID, c3InitCon);
	//spec.setBoundaryCondition("periodic","south", c4ID, c4InitCon);
	//spec.setBoundaryCondition("periodic","south", c5ID, c5InitCon);
	//spec.setBoundaryCondition("periodic","south", c6ID, c6InitCon);

	//spec.setBoundaryCondition("periodic","north", c1ID, c1InitCon);
	//spec.setBoundaryCondition("periodic","north", c2ID, c2InitCon);
	//spec.setBoundaryCondition("periodic","north", c3ID, c3InitCon);
	//spec.setBoundaryCondition("periodic","north", c4ID, c4InitCon);
	//spec.setBoundaryCondition("periodic","north", c5ID, c5InitCon);
	//spec.setBoundaryCondition("periodic","north", c6ID, c6InitCon);

	// Set source
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){

			meshCell* cell = model.getCellByLoc(i,j);
			y = cell->y;
			y1 = y - model.dy/2.;
			y2 = y + model.dy/2.;
			s = (1./model.dy)*(yLength/M_PI)*(cos(M_PI*y1/yLength) - 
				cos(M_PI*y2/yLength));

			//std::cout << y1 << " " << y2 << " " << s << std::endl;
			spec.setSpeciesSource(i, j, c1ID, c1Coeffs, a1);
			spec.setSpeciesSource(i, j, c2ID, c2Coeffs, a2);
			spec.setSpeciesSource(i, j, c3ID, c3Coeffs, a3);
			spec.setSpeciesSource(i, j, c4ID, c4Coeffs, a4);
			spec.setSpeciesSource(i, j, c5ID, c5Coeffs, a5);
			spec.setSpeciesSource(i, j, c6ID, c6Coeffs, a6);
		}
	}

	for (int k = 0; k < steps; k++){
		t = t + dt;
		// Solve with CRAM
		spec.solve(t);
		//spec.solveImplicit(t);
		//spec.solve();

		outputFile.open("precursorsSingleChan.out", std::ios_base::app);
		outputFile << "Time: "+std::to_string(t)+"\n";
		//printf (" %4.6f \n", t);
		// Gets species Concentrations
	}
	//std::cout.precision(16);
		if (myid==0){
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					meshCell* cell = model.getCellByLoc(i,j);
					y = cell->y;
					if (yVelocity != 0.0){
						c1Ana = precursorAnalitical(y, yVelocity, a1, yLength, lambdaC1);
						c2Ana = precursorAnalitical(y, yVelocity, a2, yLength, lambdaC2);
						c3Ana = precursorAnalitical(y, yVelocity, a3, yLength, lambdaC3);
						c4Ana = precursorAnalitical(y, yVelocity, a4, yLength, lambdaC4);
						c5Ana = precursorAnalitical(y, yVelocity, a5, yLength, lambdaC5);
						c6Ana = precursorAnalitical(y, yVelocity, a6, yLength, lambdaC6);
					}
					else{
						c1Ana = precursorAnalitical(t, yVelocity, a1, yLength, lambdaC1);
						c2Ana = precursorAnalitical(t, yVelocity, a2, yLength, lambdaC2);
						c3Ana = precursorAnalitical(t, yVelocity, a3, yLength, lambdaC3);
						c4Ana = precursorAnalitical(t, yVelocity, a4, yLength, lambdaC4);
						c5Ana = precursorAnalitical(t, yVelocity, a5, yLength, lambdaC5);
						c6Ana = precursorAnalitical(t, yVelocity, a6, yLength, lambdaC6);
					}

					c1Con = spec.getSpecies(i, j, c1ID);
					c2Con = spec.getSpecies(i, j, c2ID);
					c3Con = spec.getSpecies(i, j, c3ID);
					c4Con = spec.getSpecies(i, j, c4ID);
					c5Con = spec.getSpecies(i, j, c5ID);
					c6Con = spec.getSpecies(i, j, c6ID);

					c1Error = std::abs(c1Con-c1Ana);
					c2Error = std::abs(c2Con-c2Ana);
					c3Error = std::abs(c3Con-c3Ana);
					c4Error = std::abs(c4Con-c4Ana);
					c5Error = std::abs(c5Con-c5Ana);
					c6Error = std::abs(c6Con-c6Ana);

					c1Errorl1 = std::max(c1Errorl1, c1Error);
					c2Errorl1 = std::max(c2Errorl1, c2Error);
					c3Errorl1 = std::max(c3Errorl1, c3Error);
					c4Errorl1 = std::max(c4Errorl1, c4Error);
					c5Errorl1 = std::max(c5Errorl1, c5Error);
					c6Errorl1 = std::max(c6Errorl1, c6Error);

					//c1Error = 0.0;
					//c2Error = 0.0;
					//c3Error = 0.0;
					//c4Error = 0.0;
					//c5Error = std::abs(c5Con-c5Ana)/c5Ana;
					//c6Error = std::abs(c6Con-c6Ana)/c6Ana;

					//printf (" %2i %2i %e %e %e %e %e %e \n", i, j, c1Error, c2Error, c3Error, c4Error, c5Error,
					//c6Error);
					//	c2Con, c3Con, c4Con, c5Con, c6Con);
					//outputFile << i << " " << j << " " << c1Con << " " << c2Con << " " 
					//std::cout << i << " " << j << " " << c1Con << " " << c2Con << " " 
					//<< c3Con << " " << c4Con << " " << c5Con << " " << c6Con << std::endl;

				}
			}
		}
		//std::cout << " " << std::endl;
		//std::cout.precision(16);
		//printf (" %e %e %e %e %e %e \n", c1Errorl1, c2Errorl1, 
		//c3Errorl1, c4Errorl1, c5Errorl1, c6Errorl1);
		//std::cout << c1Errorl1 << " " << c2Errorl1 << " " << c3Errorl1 << " " << 
		//c4Errorl1 << " " << c5Errorl1 << " " << c6Errorl1 << std::endl;
		//spec.resetMatrix();
		//std::cout << " " << std::endl;
	//}

	
	model.clean();
	spec.clean();

}
//*****************************************************************************
// Test neutron precursors for multi chans
//*****************************************************************************
void testNeutronPrecursorsMultiChanFlow(int myid){
	double t = 0.0;
	int steps = 2;
	double totalTime = 10.0;
	double dt = totalTime/steps;
	int xCells = 14, yCells = 32;
	double xLength = 4.5, yLength = 15.5;
	double scale;
	double AvogNum = 6.02214076E23;
	int c1ID, c2ID, c3ID, c4ID, c5ID, c6ID;
	double c1Con, c2Con, c3Con, c4Con, c5Con, c6Con;
   double lambdaC1 = -0.0125, lambdaC2 = -0.0318, lambdaC3 = -0.109;
	double lambdaC4 = -0.3170, lambdaC5 = -1.3500, lambdaC6 = -8.640;
	double c1InitCon = 0.0, c2InitCon = 0.0, c3InitCon = 0.0, c4InitCon = 0.0;
	double c5InitCon = 0.0, c6InitCon = 0.0;
	double D_c1 = 0.0, D_c2 = 0.0, D_c3 = 0.0, D_c4 = 0.0, D_c5 = 0.0, D_c6 = 0.0;
	double s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, s5 = 0.0, s6 = 0.0;
	double maxC1 = 0.0, maxC2 = 0.0, maxC3 = 0.0, maxC4 = 0.0, maxC5 = 0.0, maxC6 = 0.0;
	double velocityScale = 1.0;
   MatrixD coeff(16,7);
	Tensor<double, 3> ceoff3d(6, 13, 16);
	std::ofstream outputFile;
	std::vector<double> c1Coeffs = {lambdaC1, 0.0, 0.0, 0.0, 0.0, 0.0};
	std::vector<double> c2Coeffs = {0.0, lambdaC2, 0.0, 0.0, 0.0, 0.0};
	std::vector<double> c3Coeffs = {0.0, 0.0, lambdaC3, 0.0, 0.0, 0.0};
	std::vector<double> c4Coeffs = {0.0, 0.0, 0.0, lambdaC4, 0.0, 0.0};
	std::vector<double> c5Coeffs = {0.0, 0.0, 0.0, 0.0, lambdaC5, 0.0};
	std::vector<double> c6Coeffs = {0.0, 0.0, 0.0, 0.0, 0.0, lambdaC6};

   // Coefficients in order of precursor groups for columns
	coeff(0,0)  = 1.9490E-04, coeff(0,1)  = 1.0149E-03, coeff(0,2)  = 9.8409E-04;
	coeff(1,0)  = 3.5140E-04, coeff(1,1)  = 1.8298E-03, coeff(1,2)  = 1.7743E-03;
	coeff(2,0)  = 4.8659E-04, coeff(2,1)  = 2.5338E-03, coeff(2,2)  = 2.4569E-03;
	coeff(3,0)  = 6.0108E-04, coeff(3,1)  = 3.1300E-03, coeff(3,2)  = 3.0350E-03;
	coeff(4,0)  = 6.9421E-04, coeff(4,1)  = 3.6150E-03, coeff(4,2)  = 3.5052E-03;
	coeff(5,0)  = 7.6455E-04, coeff(5,1)  = 3.9812E-03, coeff(5,2)  = 3.8604E-03;
	coeff(6,0)  = 8.1053E-04, coeff(6,1)  = 4.2207E-03, coeff(6,2)  = 4.0926E-03;
	coeff(7,0)  = 8.3098E-04, coeff(7,1)  = 4.3271E-03, coeff(7,2)  = 4.1958E-03;
	coeff(8,0)  = 8.2533E-04, coeff(8,1)  = 4.2977E-03, coeff(8,2)  = 4.1673E-03;
	coeff(9,0)  = 7.9376E-04, coeff(9,1)  = 4.1333E-03, coeff(9,2)  = 4.0079E-03;
	coeff(10,0) = 7.3720E-04, coeff(10,1) = 3.8388E-03, coeff(10,2) = 3.7223E-03;
	coeff(11,0) = 6.5731E-04, coeff(11,1) = 3.4228E-03, coeff(11,2) = 3.3189E-03;
	coeff(12,0) = 5.5629E-04, coeff(12,1) = 2.8967E-03, coeff(12,2) = 2.8088E-03;
	coeff(13,0) = 4.3661E-04, coeff(13,1) = 2.2735E-03, coeff(13,2) = 2.2045E-03;
	coeff(14,0) = 3.0078E-04, coeff(14,1) = 1.5662E-03, coeff(14,2) = 1.5187E-03;
	coeff(15,0) = 1.5123E-04, coeff(15,1) = 7.8752E-04, coeff(15,2) = 7.6362E-04;

   coeff(0,3)  = 2.8054E-03, coeff(0,4)  = 8.1466E-04, coeff(0,5)  = 2.8792E-04;
	coeff(1,3)  = 5.0582E-03, coeff(1,4)  = 1.4688E-03, coeff(1,5)  = 5.1911E-04;
	coeff(2,3)  = 7.0042E-03, coeff(2,4)  = 2.0339E-03, coeff(2,5)  = 7.1883E-04;
	coeff(3,3)  = 8.6521E-03, coeff(3,4)  = 2.5124E-03, coeff(3,5)  = 8.8795E-04;
	coeff(4,3)  = 9.9927E-03, coeff(4,4)  = 2.9017E-03, coeff(4,5)  = 1.0255E-03;
	coeff(5,3)  = 1.1005E-02, coeff(5,4)  = 3.1957E-03, coeff(5,5)  = 1.1294E-03;
	coeff(6,3)  = 1.1667E-02, coeff(6,4)  = 3.3880E-03, coeff(6,5)  = 1.1974E-03;
	coeff(7,3)  = 1.1961E-02, coeff(7,4)  = 3.4734E-03, coeff(7,5)  = 1.2276E-03;
	coeff(8,3)  = 1.1880E-02, coeff(8,4)  = 3.4498E-03, coeff(8,5)  = 1.2192E-03;
	coeff(9,3)  = 1.1426E-02, coeff(9,4)  = 3.3178E-03, coeff(9,5)  = 1.1726E-03;
	coeff(10,3) = 1.0611E-02, coeff(10,4) = 3.0814E-03, coeff(10,5) = 1.0890E-03;
	coeff(11,3) = 9.4616E-03, coeff(11,4) = 2.7475E-03, coeff(11,5) = 9.7103E-04;
	coeff(12,3) = 8.0074E-03, coeff(12,4) = 2.3252E-03, coeff(12,5) = 8.2178E-04;
	coeff(13,3) = 6.2846E-03, coeff(13,4) = 1.8250E-03, coeff(13,5) = 6.4498E-04;
	coeff(14,3) = 4.3295E-03, coeff(14,4) = 1.2572E-03, coeff(14,5) = 4.4433E-04;
	coeff(15,3) = 2.1769E-03, coeff(15,4) = 6.3215E-04, coeff(15,5) = 2.2341E-04;


	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);

	// Sets the y velocity
	model.setConstantYVelocity(velocityScale*2.16*3.28, 0);
	model.setConstantYVelocity(velocityScale*2.34*3.38, 1);
	model.setConstantYVelocity(velocityScale*2.46*3.38, 2);
	model.setConstantYVelocity(velocityScale*2.24*3.38, 3);
	model.setConstantYVelocity(velocityScale*2.1*3.38, 4);
	model.setConstantYVelocity(velocityScale*2.05*3.38, 5);
	model.setConstantYVelocity(velocityScale*1.95*3.38, 6);
	model.setConstantYVelocity(velocityScale*1.93*3.38, 7);
	model.setConstantYVelocity(velocityScale*1.9*3.38, 8);
	model.setConstantYVelocity(velocityScale*1.88*3.38, 9);
	model.setConstantYVelocity(velocityScale*1.87*3.38, 10);
	model.setConstantYVelocity(velocityScale*1.86*3.38, 11);
	model.setConstantYVelocity(velocityScale*1.85*3.38, 12);
	model.setConstantYVelocity(velocityScale*1.84*3.38, 13);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds species
	c1ID = spec.addSpecies(1.0, 0.0, D_c1);
	c2ID = spec.addSpecies(1.0, 0.0, D_c2);
	c3ID = spec.addSpecies(1.0, 0.0, D_c3);
	c4ID = spec.addSpecies(1.0, 0.0, D_c4);
	c5ID = spec.addSpecies(1.0, 0.0, D_c5);
	c6ID = spec.addSpecies(1.0, 0.0, D_c6);

	// Sets BCs
	spec.setBoundaryCondition("periodic","south", c1ID, c1InitCon);
	spec.setBoundaryCondition("periodic","south", c2ID, c2InitCon);
	spec.setBoundaryCondition("periodic","south", c3ID, c3InitCon);
	spec.setBoundaryCondition("periodic","south", c4ID, c4InitCon);
	spec.setBoundaryCondition("periodic","south", c5ID, c5InitCon);
	spec.setBoundaryCondition("periodic","south", c6ID, c6InitCon);

	spec.setBoundaryCondition("periodic","north", c1ID, c1InitCon);
	spec.setBoundaryCondition("periodic","north", c2ID, c2InitCon);
	spec.setBoundaryCondition("periodic","north", c3ID, c3InitCon);
	spec.setBoundaryCondition("periodic","north", c4ID, c4InitCon);
	spec.setBoundaryCondition("periodic","north", c5ID, c5InitCon);
	spec.setBoundaryCondition("periodic","north", c6ID, c6InitCon);

	// Set source
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			meshCell* cell = model.getCellByLoc(i,j);
			scale = cos(cell->x/3.5);
			if (j < 16){
				// need to convert g/s/cm^3 to lbm/s/ft^3 by 62.427961
				s1 = scale*coeff(j,c1ID); 
				s2 = scale*coeff(j,c2ID);
				s3 = scale*coeff(j,c3ID);
				s4 = scale*coeff(j,c4ID);
				s5 = scale*coeff(j,c5ID);
				s6 = scale*coeff(j,c6ID);
			}
			else{
				s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, s5 = 0.0, s6 = 0.0;
			}
			spec.setSpeciesSource(i, j, c1ID, c1Coeffs, s1);
			spec.setSpeciesSource(i, j, c2ID, c2Coeffs, s2);
			spec.setSpeciesSource(i, j, c3ID, c3Coeffs, s3);
			spec.setSpeciesSource(i, j, c4ID, c4Coeffs, s4);
			spec.setSpeciesSource(i, j, c5ID, c5Coeffs, s5);
			spec.setSpeciesSource(i, j, c6ID, c6Coeffs, s6);
		}
	}

	for (int k = 0; k < steps; k++){
		t = t + dt;
		// Solve with CRAM
		//spec.solveImplicit(t);
		spec.solve(t);
		//spec.solve();

		//outputFile.open("precursorsMultiChan.out", std::ios_base::app);
		//outputFile << "Time: "+std::to_string(t)+"\n";
		printf (" %4.6f \n", t);
	}
	std::cout.precision(16);
		// Gets species Concentrations
		if (myid==0){
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					c1Con = spec.getSpecies(i, j, c1ID);
					c2Con = spec.getSpecies(i, j, c2ID);
					c3Con = spec.getSpecies(i, j, c3ID);
					c4Con = spec.getSpecies(i, j, c4ID);
					c5Con = spec.getSpecies(i, j, c5ID);
					c6Con = spec.getSpecies(i, j, c6ID);

					maxC1 = std::max(c1Con, maxC1);
					maxC2 = std::max(c2Con, maxC2);
					maxC3 = std::max(c3Con, maxC3);
					maxC4 = std::max(c4Con, maxC4);
					maxC5 = std::max(c5Con, maxC5);
					maxC6 = std::max(c6Con, maxC6);

					//printf (" %2i %2i %*E %*E %*E %*E %*E %*E\n", i, j, c1Con, 
					//	c2Con, c3Con, c4Con, c5Con, c6Con);
					//outputFile << i << " " << j << " " << c1Con << " " << c2Con << " " 
					std::cout << i << " " << j << " " << c1Con << " " << c2Con << " " 
						<< c3Con << " " << c4Con << " " << c5Con << " " << c6Con << std::endl;

				}
			}
		}
	//}
	//std::cout << maxC1 << std::endl;
	//std::cout << maxC2 << std::endl;
	//std::cout << maxC3 << std::endl;
	//std::cout << maxC4 << std::endl;
	//std::cout << maxC5 << std::endl;
	//std::cout << maxC6 << std::endl;

	
	model.clean();
	spec.clean();

}
//*****************************************************************************
// Test benchmark for Ben
//*****************************************************************************
void testBenBenchmark(int myid){
	int xCells = 1, yCells = 100;
	double xLength = 1.0, yLength = 300.;
	double yVelocity = 30.0;
	double specInitCon = 0.0;
	double specMM = 1.0;
	double t = 10000000.0;
   double lambda_spec = -.3;
	double D_spec = 1.0;
	double specS = 0.0, x = 0.0;
	meshCell* cell;
	int specID, specID2;
	double specCon;
	std::vector<double> specCoeffs = {lambda_spec};

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);

	// Sets the x velocity
	model.setConstantYVelocity(yVelocity);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds xenon and iodine species
	specID = spec.addSpecies(specMM, 0.0, D_spec);
	spec.setBoundaryCondition("dirichlet", "south", specID, 1.0);

	// Set source
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			cell = model.getCellByLoc(i,j);
			x = cell->x;
			specS = (x < 100.) ? sin(M_PI*x/100.) : 0.0;
			//std::cout << x << " " << specS << std::endl;
			spec.setSpeciesSource(i, j, specID, specCoeffs, specS);
		}
	}

	// Solve with CRAM
	spec.solve(t);

	std::ofstream outputFile;
	outputFile.open("benBenchmark.out", std::ios_base::app);
	outputFile << "Time: "+std::to_string(t)+"\n";
	// Gets species Concentrations
	if (myid==0){
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				specCon = spec.getSpecies(i, j, specID);

				//std::cout << specCon << std::endl;
				//printf (" %2i %2i %E \n", i, j, specCon);
				outputFile << i << " " << j << " " << specCon << std::endl;
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

	testXenonIodineNoFlow(myid);
	testProblem1NoFlow(myid);
	testXenonIodineYFlow(myid);
	testXenonIodineXFlow(myid);
	testDiffusion2D(myid);
	testNeutronPrecursorsFlow(myid);
	testNeutronPrecursorsMultiChanFlow(myid);
	testBenBenchmark(myid);

	mpi.finalize();
}
