#include <Eigen/Core>
#include <Eigen/Sparse>
#include <assert.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "CRAM.h"
#include "mpiProcess.h"
#include "modelMesh.h"

using namespace Eigen;

void testInit(){

	modelMesh mesh(2, 5, 1.0, 1.0);
	mesh.setConstantXVelocity(2.0, 0);
	mesh.setConstantXVelocity(2.0, 1);
	mesh.setConstantXVelocity(2.0);
	mesh.setConstantYVelocity(0.1);
	mesh.setConstantYVelocity(0.2, 4);
	mesh.clean();

}
void testInitSpecies(){
	int xCells = 2, yCells = 5;
	double xLength = 1.0, yLength = 1.0;
	double spec1InitCon = 1.0, spec2InitCon = 2.0;
	double spec1MM = 2.0, spec2MM = 3.0;
	int specID1, specID2;
	double spec1Con, spec2Con;

	modelMesh model(xCells, yCells, xLength, yLength);
	specID1 = model.addSpecies(spec1MM, spec1InitCon);
	specID2 = model.addSpecies(spec2MM, spec2InitCon);

	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			// Gets species pointers
			species* spec1 = model.getSpeciesPtr(i, j, specID1);
			species* spec2 = model.getSpeciesPtr(i, j, specID2);

			// Gets species Concentrations
			spec1Con = model.getSpecies(i, j, specID1);;
			spec2Con = model.getSpecies(i, j, specID2);;

			// Makes sure all species concentrations are right
			assert(1.0 == spec1Con);
			assert(2.0 == spec2Con);
			// Makes sure all species molar masses are right
			assert(spec1->MM == spec1MM);
			assert(spec2->MM == spec2MM);
		}
	}
	
	model.clean();

}

int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	testInit();
	testInitSpecies();
}
