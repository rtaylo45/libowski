#include <Eigen/Core>
#include <Eigen/Sparse>
#include "coreSolver.h"
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

	modelMesh dude(1, 5, 1.0, 1.0);
	dude.setConstantXVelocity(2.0);
}
