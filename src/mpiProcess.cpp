#include <iostream>
#include "mpiProcess.h"

// Global mpi object
mpiProcess mpi;

void mpiProcess::initMPI(){
#ifdef HAVE_MPI
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	isInit = true;
#endif
}

void mpiProcess::finalize(){
#ifdef HAVE_MPI
	MPI_Finalize();
#endif
}

void mpiProcess::send(VectorXcd x, int count, int id, int MTAG){
#ifdef HAVE_MPI
	MPI_Send(x.data(), count, MPI::DOUBLE_COMPLEX, id, MTAG, MPI_COMM_WORLD);
#endif

}
//**************************************************************************
// Receives  MPI data
//**************************************************************************
VectorXcd mpiProcess::recv(VectorXcd x, int count, int islave, int MTAG){
	VectorXcd xRet;
#ifdef HAVE_MPI
	MPI_Status status;
	MPI_Recv(x.data(), count, MPI::DOUBLE_COMPLEX, islave, MTAG, 
		MPI_COMM_WORLD, &status);
	xRet = x;
#endif
	return xRet;

}

