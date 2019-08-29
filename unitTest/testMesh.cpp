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
	int specID1, specID2;
	modelMesh mesh(2, 5, 1.0, 1.0);
	specID1 = mesh.addSpecies(1.0, 2.0);
	specID2 = mesh.addSpecies(2.0, 2.0);
	mesh.clean();

}

int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	testInit();
	testInitSpecies();
}
