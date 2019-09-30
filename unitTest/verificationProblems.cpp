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

#include "CRAM.h"
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

	double diff = abs(goalVal - testVal);
	if (diff < rtol) { retBool = true; }
	if (diff/goalVal < atol) { retBool = true; }
	return retBool;
}
//*****************************************************************************
// Test that the species driver sets up the problem right and solves the 
// system right. 
//*****************************************************************************
void testXenonIodineNoFlow(){
	int xCells = 1, yCells = 10;
	double xLength = 1.0, yLength = 1.0;
	double xenonInitCon = 0.0, iodineInitCon = 0.0;
	double xenonMM = 135.0, iodineMM = 135.0;
	double numOfSteps = 10.;
	double tEnd = 10000.0;
	double t;
	double dt = tEnd/numOfSteps;
   double lambda_I = 2.11E-5;
   double lambda_xe = 2.9306E-5;
   double sigma_a = 2.002E-22;
   double Sigma_f = 9.7532E-1;
   double flux = 2.5E16;
   double gamma_xe = 0.002468;
   double gamma_I = 0.063033;
	double N_xe_0 = 0.0, N_I_0 = 0.0;
	int xenonID, iodineID;
	double xenonCon, iodineCon;
	std::vector<double> xenonCoeffs = {-lambda_xe-sigma_a*flux, lambda_I};
	std::vector<double> iodineCoeffs = {0.0, -lambda_I};
	double xenonS = gamma_xe*Sigma_f*flux;
	double iodineS = gamma_I*Sigma_f*flux;

	modelMesh model(xCells, yCells, xLength, yLength);
	speciesDriver spec = speciesDriver(&model);
	xenonID = spec.addSpecies(xenonMM, N_xe_0);
	iodineID = spec.addSpecies(iodineMM, N_I_0);

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
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				xenonCon = spec.getSpecies(i, j, xenonID);
				iodineCon = spec.getSpecies(i, j, iodineID);
				assert(isApprox(xenonCon, N_xe));
				assert(isApprox(iodineCon, N_I));
			}
		}
	}
	
	model.clean();
	spec.clean();
}
//*****************************************************************************
// Test Xenon iodine flow problem
//*****************************************************************************
void testXenonIodineFlow(){
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
	xenonID = spec.addSpecies(xenonMM, N_xe_0);
	iodineID = spec.addSpecies(iodineMM, N_I_0);
	spec.setBoundaryCondition(0, 0, xenonID, xenonInitCon);
	spec.setBoundaryCondition(0, 0, iodineID, iodineInitCon);

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
			error = std::max(std::abs(iodineCon - N_I)/N_I, error);
			//std::cout << y << " " << iodineCon << " " << N_I << " " << std::abs(iodineCon - N_I)/N_I << std::endl;
			assert(isApprox(iodineCon, N_I));
		}
	}
	//std::cout << "Max l-1 error: " << error << std::endl;

	
	model.clean();
	spec.clean();

}
//*****************************************************************************
// Test neutron precursors for sinlge channel
//*****************************************************************************
void testNeutronPrecursorsFlow(){
	double t = 0.0;
	int steps = 2;
	double totalTime = 140.0;
	double dt = totalTime/steps;
	int xCells = 1, yCells = 16;
	double xLength = 1.0, yLength = 8.0;
	double scale;
	double AvogNum = 6.02214076E23;
	int c1ID, c2ID, c3ID, c4ID, c5ID, c6ID;
	double c1Con, c2Con, c3Con, c4Con, c5Con, c6Con;
   double lambdaC1 = -0.0125, lambdaC2 = -0.0318, lambdaC3 = -0.109;
	double lambdaC4 = -0.3170, lambdaC5 = -1.3500, lambdaC6 = -8.640;
	double c1InitCon = 0.0, c2InitCon = 0.0, c3InitCon = 0.0, c4InitCon = 0.0;
	double c5InitCon = 0.0, c6InitCon = 0.0;
   MatrixXd coeff(16,7);
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
	model.setConstantYVelocity(1.0, 0);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds species
	c1ID = spec.addSpecies(1.0, 0.0);
	c2ID = spec.addSpecies(1.0, 0.0);
	c3ID = spec.addSpecies(1.0, 0.0);
	c4ID = spec.addSpecies(1.0, 0.0);
	c5ID = spec.addSpecies(1.0, 0.0);
	c6ID = spec.addSpecies(1.0, 0.0);

	// Sets BCs
	for (int i = 0; i < xCells; i++){
		spec.setBoundaryCondition(i, 0, c1ID, c1InitCon);
		spec.setBoundaryCondition(i, 0, c2ID, c2InitCon);
		spec.setBoundaryCondition(i, 0, c3ID, c3InitCon);
		spec.setBoundaryCondition(i, 0, c4ID, c4InitCon);
		spec.setBoundaryCondition(i, 0, c5ID, c5InitCon);
		spec.setBoundaryCondition(i, 0, c6ID, c6InitCon);
	}

	// Set source
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			spec.setSpeciesSource(i, j, c1ID, c1Coeffs, coeff(j,c1ID));
			spec.setSpeciesSource(i, j, c2ID, c2Coeffs, coeff(j,c2ID));
			spec.setSpeciesSource(i, j, c3ID, c3Coeffs, coeff(j,c3ID));
			spec.setSpeciesSource(i, j, c4ID, c4Coeffs, coeff(j,c4ID));
			spec.setSpeciesSource(i, j, c5ID, c5Coeffs, coeff(j,c5ID));
			spec.setSpeciesSource(i, j, c6ID, c6Coeffs, coeff(j,c6ID));
		}
	}

	for (int k = 0; k < steps; k++){
		t = t + dt;
		// Solve with CRAM
		spec.solve(t);

		std::ofstream outputFile;
		outputFile.open("precursorsSingleChan.out", std::ios_base::app);
		outputFile << "Time: "+std::to_string(t)+"\n";
		//printf (" %4.6f \n", t);
		// Gets species Concentrations
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				c1Con = spec.getSpecies(i, j, c1ID);
				c2Con = spec.getSpecies(i, j, c2ID);
				c3Con = spec.getSpecies(i, j, c3ID);
				c4Con = spec.getSpecies(i, j, c4ID);
				c5Con = spec.getSpecies(i, j, c5ID);
				c6Con = spec.getSpecies(i, j, c6ID);

				//printf (" %2i %2i %E %E %E %E %E %E\n", i, j, c1Con, 
				//	c2Con, c3Con, c4Con, c5Con, c6Con);
				outputFile << i << " " << j << " " << c1Con << " " << c2Con << " " 
					<< c3Con << " " << c4Con << " " << c5Con << " " << c6Con << std::endl;

			}
		}
		spec.resetMatrix();
		//std::cout << " " << std::endl;
	}

	
	model.clean();
	spec.clean();

}
//*****************************************************************************
// Test neutron precursors for multi chans
//*****************************************************************************
void testNeutronPrecursorsMultiChanFlow(){
	double t = 0.0;
	int steps = 100;
	double totalTime = 5.0;
	double dt = totalTime/steps;
	int xCells = 14, yCells = 16;
	double xLength = 4.5, yLength = 5.5;
	double scale;
	double AvogNum = 6.02214076E23;
	int c1ID, c2ID, c3ID, c4ID, c5ID, c6ID;
	double c1Con, c2Con, c3Con, c4Con, c5Con, c6Con;
   double lambdaC1 = -0.0125, lambdaC2 = -0.0318, lambdaC3 = -0.109;
	double lambdaC4 = -0.3170, lambdaC5 = -1.3500, lambdaC6 = -8.640;
	double c1InitCon = 0.0, c2InitCon = 0.0, c3InitCon = 0.0, c4InitCon = 0.0;
	double c5InitCon = 0.0, c6InitCon = 0.0;
   MatrixXd coeff(16,7);
	Tensor<double, 3> ceoff3d(6, 13, 16);
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
	model.setConstantYVelocity(2.16*3.28, 0);
	model.setConstantYVelocity(2.34*3.38, 1);
	model.setConstantYVelocity(2.46*3.38, 2);
	model.setConstantYVelocity(2.24*3.38, 3);
	model.setConstantYVelocity(2.1*3.38, 4);
	model.setConstantYVelocity(2.05*3.38, 5);
	model.setConstantYVelocity(1.95*3.38, 6);
	model.setConstantYVelocity(1.93*3.38, 7);
	model.setConstantYVelocity(1.9*3.38, 8);
	model.setConstantYVelocity(1.88*3.38, 9);
	model.setConstantYVelocity(1.87*3.38, 10);
	model.setConstantYVelocity(1.86*3.38, 11);
	model.setConstantYVelocity(1.85*3.38, 12);
	model.setConstantYVelocity(1.84*3.38, 13);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds species
	c1ID = spec.addSpecies(1.0, 0.0);
	c2ID = spec.addSpecies(1.0, 0.0);
	c3ID = spec.addSpecies(1.0, 0.0);
	c4ID = spec.addSpecies(1.0, 0.0);
	c5ID = spec.addSpecies(1.0, 0.0);
	c6ID = spec.addSpecies(1.0, 0.0);

	// Sets BCs
	for (int i = 0; i < xCells; i++){
		spec.setBoundaryCondition(i, 0, c1ID, c1InitCon);
		spec.setBoundaryCondition(i, 0, c2ID, c2InitCon);
		spec.setBoundaryCondition(i, 0, c3ID, c3InitCon);
		spec.setBoundaryCondition(i, 0, c4ID, c4InitCon);
		spec.setBoundaryCondition(i, 0, c5ID, c5InitCon);
		spec.setBoundaryCondition(i, 0, c6ID, c6InitCon);
	}

	// Set source
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			meshCell* cell = model.getCellByLoc(i,j);
			scale = cos(cell->x/3.5);
			spec.setSpeciesSource(i, j, c1ID, c1Coeffs, scale*coeff(j,c1ID));
			spec.setSpeciesSource(i, j, c2ID, c2Coeffs, scale*coeff(j,c2ID));
			spec.setSpeciesSource(i, j, c3ID, c3Coeffs, scale*coeff(j,c3ID));
			spec.setSpeciesSource(i, j, c4ID, c4Coeffs, scale*coeff(j,c4ID));
			spec.setSpeciesSource(i, j, c5ID, c5Coeffs, scale*coeff(j,c5ID));
			spec.setSpeciesSource(i, j, c6ID, c6Coeffs, scale*coeff(j,c6ID));
		}
	}

	for (int k = 0; k < steps; k++){
		t = t + dt;
		// Solve with CRAM
		spec.solve(t);

		std::ofstream outputFile;
		outputFile.open("precursorsMultiChan.out", std::ios_base::app);
		outputFile << "Time: "+std::to_string(t)+"\n";
		printf (" %4.6f \n", t);
		// Gets species Concentrations
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				c1Con = spec.getSpecies(i, j, c1ID);
				c2Con = spec.getSpecies(i, j, c2ID);
				c3Con = spec.getSpecies(i, j, c3ID);
				c4Con = spec.getSpecies(i, j, c4ID);
				c5Con = spec.getSpecies(i, j, c5ID);
				c6Con = spec.getSpecies(i, j, c6ID);

				//printf (" %2i %2i %E %E %E %E %E %E\n", i, j, c1Con, 
				//	c2Con, c3Con, c4Con, c5Con, c6Con);
				outputFile << i << " " << j << " " << c1Con << " " << c2Con << " " 
					<< c3Con << " " << c4Con << " " << c5Con << " " << c6Con << std::endl;

			}
		}
		//std::cout << " " << std::endl;
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

	testXenonIodineNoFlow();
	testXenonIodineFlow();
	testNeutronPrecursorsFlow();
	testNeutronPrecursorsMultiChanFlow();

	mpi.finalize();
}
