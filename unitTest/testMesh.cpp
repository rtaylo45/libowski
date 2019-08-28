#include <Eigen/Core>
#include <Eigen/Sparse>
#include "CRAM.h"
#include "mpiProcess.h"
#include "modelMesh.h"
#include <assert.h>
#include <iostream>
#include <vector>
#include <random>
#include <math.h>

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

int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	testInit();
}
