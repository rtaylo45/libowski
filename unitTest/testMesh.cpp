#include <Eigen/Core>
#include <Eigen/Sparse>
#include <assert.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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
// test the init of the mesh class. For now i kinda just test to make sure
// the function run.
//*****************************************************************************
void testInit(){

	modelMesh mesh(2, 5, 1.0, 1.0);
	mesh.setConstantXVelocity(2.0, 0);
	mesh.setConstantXVelocity(2.0, 1);
	mesh.setConstantXVelocity(2.0);
	mesh.setConstantYVelocity(0.1);
	mesh.setConstantYVelocity(0.2, 1);
	mesh.clean();

}

//*****************************************************************************
// test the init of the species driver. Test to make sure multiple species
// can be added. initial concentrations are added correctly. Molar masses
// are added currect and sources are set right. 
//*****************************************************************************
void testSpeciesDriver(){
	int xCells = 2, yCells = 5;
	double xLength = 1.0, yLength = 1.0;
	double spec1InitCon = 1.0, spec2InitCon = 2.0;
	double spec1MM = 2.0, spec2MM = 3.0;
	int specID1, specID2;
	double spec1Con, spec2Con;
	std::vector<double> spec1Coeffs = {5.0, 10.0};
	std::vector<double> spec2Coeffs = {6.0, 11.0};
	double spec1S = 50.0, spec2S = 100.0;

	modelMesh model(xCells, yCells, xLength, yLength);
	speciesDriver spec = speciesDriver(&model);
	specID1 = spec.addSpecies(spec1MM, spec1InitCon);
	specID2 = spec.addSpecies(spec2MM, spec2InitCon);

	// Sets sourse terms for the model
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			spec.setSpeciesSource(i,j,specID1,	spec1Coeffs, spec1S);
			spec.setSpeciesSource(i,j,specID2,	spec2Coeffs, spec2S);
		}
	}
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			// Gets species pointers
			species* spec1 = spec.getSpeciesPtr(i, j, specID1);
			species* spec2 = spec.getSpeciesPtr(i, j, specID2);

			// Gets species Concentrations
			spec1Con = spec.getSpecies(i, j, specID1);;
			spec2Con = spec.getSpecies(i, j, specID2);;

			// Makes sure all species concentrations are right
			assert(1.0 == spec1Con);
			assert(2.0 == spec2Con);
			// Makes sure all species molar masses are right
			assert(spec1->MM == spec1MM);
			assert(spec2->MM == spec2MM);
		}
	}
	
	model.clean();
	spec.clean();

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
// Test that the species driver sets up the problem right and solves the 
// system right. 
//*****************************************************************************
void testYDirectionAdvection(){
	int xCells = 1, yCells = 10;
	double xLength = 1.0, yLength = 100.0;
	double specInitCon = 10.0;
	double numOfSteps = 5.;
	double tEnd = 60.0; 
	double dt = tEnd/numOfSteps;
	double yVelocity = 6.0; // ft/s
	int specID, spec2ID;
	std::vector<double> coeffs = {0.0};
	double dx = yLength/(float)yCells;

	// Builds the mesh
	modelMesh model(xCells, yCells, xLength, yLength);

	// Sets the x velocity
	model.setConstantYVelocity(yVelocity);

	// Sets species driver
	speciesDriver spec = speciesDriver(&model);

	// Adds species to the model
	specID = spec.addSpecies(1.0, 10.0);
	spec.setBoundaryCondition(0, 0, specID, 2.*specInitCon);

	for (int step = 1; step <= numOfSteps; step++){
		double t = step*dt;
		// Solve with CRAM
		spec.solve(t);

		// Outlet concentration	
		double specConOut = spec.getSpecies(0, yCells-1, specID);
		double specConIn = spec.getSpecies(0, 0, specID);
		printf (" %4.2f %6.3f %6.3f\n", t, specConIn, specConOut);
	}
	model.clean();
	spec.clean();
}
//*****************************************************************************
// Test Xenon iodine flow problem
//*****************************************************************************
void testXenonIodineFlow(){
	int xCells = 1, yCells = 10;
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
   //double flux = 2.5E10;
   double gamma_xe = 0.002468;
   double gamma_I = 0.063033;
	double N_xe_0 = 0.0, N_I_0 = 0.0;
	int xenonID, iodineID;
	double xenonCon, iodineCon;
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
			std::cout << y << " " << iodineCon << " " << N_I << " " << std::abs(iodineCon - N_I)/N_I << std::endl;
			//assert(isApprox(iodineCon, N_I, 1e-4, 1e-5));
		}
	}

	
	model.clean();
	spec.clean();

}
//*****************************************************************************
// Test neutron precursors
//*****************************************************************************
//void testNeutronPrecursorsFlow(){
//	int xCells = 7, yCells = 16;
//	double xLength = 2.5, yLength = 8.0;
//	double AvogNum = 6.02214076E23;
//	double t = 10000000.0;
//	int c1ID, c2ID, c3ID, c4ID, c5ID, c6ID;
//	double c1Con, c2Con, c3Con, c4Con, c5Con;
//	std::vector<double> xenonCoeffs = {-lambda_xe-sigma_a*flux, lambda_I};
//	std::vector<double> iodineCoeffs = {0.0, -lambda_I};
//	double xenonS = gamma_xe*Sigma_f*flux*xenonMM/AvogNum;
//	double iodineS = gamma_I*Sigma_f*flux*iodineMM/AvogNum;
//
//	// Builds the mesh
//	modelMesh model(xCells, yCells, xLength, yLength);
//
//	// Sets the x velocity
//	model.setConstantYVelocity(yVelocity);
//
//	// Sets species driver
//	speciesDriver spec = speciesDriver(&model);
//
//	// Adds xenon and iodine species
//	xenonID = spec.addSpecies(xenonMM, N_xe_0);
//	iodineID = spec.addSpecies(iodineMM, N_I_0);
//	spec.setBoundaryCondition(0, 0, xenonID, xenonInitCon);
//	spec.setBoundaryCondition(0, 0, iodineID, iodineInitCon);
//
//	// Set source
//	for (int i = 0; i < xCells; i++){
//		for (int j = 0; j < yCells; j++){
//			spec.setSpeciesSource(i, j, xenonID, xenonCoeffs, xenonS);
//			spec.setSpeciesSource(i, j, iodineID, iodineCoeffs, iodineS);
//		}
//	}
//
//	// Solve with CRAM
//	spec.solve(t);
//
//	// Gets species Concentrations
//	for (int i = 0; i < xCells; i++){
//		for (int j = 0; j < yCells; j++){
//			xenonCon = spec.getSpecies(i, j, xenonID);
//			iodineCon = spec.getSpecies(i, j, iodineID);
//
//			//std::cout << xenonCon << " " << iodineCon << " " << N_I << std::endl;
//			std::cout << y << " " << iodineCon << " " << N_I << " " << std::abs(iodineCon - N_I)/N_I << std::endl;
//		}
//	}
//
//	
//	model.clean();
//	spec.clean();
//
//}
//*****************************************************************************
// Main test
//*****************************************************************************
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	//testInit();
	//testSpeciesDriver();
	//testXenonIodineNoFlow();
	//testYDirectionAdvection();
	testXenonIodineFlow();
	//testNeutronPrecursorsFlow();

	mpi.finalize();
}
