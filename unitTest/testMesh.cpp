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

int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	modelMesh dude(2, 5, 1.0, 1.0);
	dude.setConstantXVelocity(2.0, 0);
	dude.setConstantXVelocity(2.0, 1);
	dude.setConstantXVelocity(2.0);
	dude.setConstantYVelocity(0.1);
	dude.setConstantYVelocity(0.2, 4);
}
