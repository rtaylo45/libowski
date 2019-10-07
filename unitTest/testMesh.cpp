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
	specID1 = spec.addSpecies(spec1MM, spec1InitCon, 0.0);
	specID2 = spec.addSpecies(spec2MM, spec2InitCon, 0.0);

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
	specID = spec.addSpecies(1.0, 10.0, 0.0);
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
// Main test
//*****************************************************************************
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	testInit();
	testSpeciesDriver();
	//testYDirectionAdvection();

	mpi.finalize();
}
