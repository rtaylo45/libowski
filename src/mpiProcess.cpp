#include "mpiProcess.h"

// Global mpi object
mpiProcess mpi;

//**************************************************************************
// Constructor
//**************************************************************************
mpiProcess::mpiProcess(){initMPI();};

//**************************************************************************
// Init the mpi process
//**************************************************************************
void mpiProcess::initMPI(){
#ifdef HAVE_MPI
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	isInit = true;
#endif
}

//**************************************************************************
// Finalize the mpi process
//**************************************************************************
void mpiProcess::finalize(){
#ifdef HAVE_MPI
	MPI_Finalize();
	isInit = false;
#endif
}

//**************************************************************************
// Sends complex long double eigen sparse matrix data
//
// @param x			Complex eigen sparse matrix to send
// @param count	Number of matrix elements
// @param id		Processor ID
// @param MTAG		Message tag
//**************************************************************************
void mpiProcess::send(const SparseMatrixCLD& x, int count, int id, int MTAG){
#ifdef HAVE_MPI
	MPI_Send(x.data(), count, MPI::LONG_DOUBLE_COMPLEX, id, MTAG, MPI_COMM_WORLD);
#endif
}

//**************************************************************************
// Receives complex long double eigen sparse matrix data
//
// @param x			Complex eigen sparse matrix to send
// @param count	Number of matrix elements
// @param id		Processor ID
// @param MTAG		Message tag
//
// MPI requires me to pass in the actual vector to recieve the data. Thats
// why the data that you recieve is an input argument to the method. If 
// HAVE_MPI is not defined then the method just returns an empty vector.
//**************************************************************************
SparseMatrixCLD mpiProcess::recv(const SparseMatrixCLD& x, int count, int islave, int MTAG){
	// Inits a zero return vector
	SparseMatrixCLD xRet;
#ifdef HAVE_MPI
	MPI_Status status;
	MPI_Recv(x.data(), count, MPI::LONG_DOUBLE_COMPLEX, islave, MTAG, 
		MPI_COMM_WORLD, &status);
	xRet = x;
#endif
	return xRet;
}

//**************************************************************************
// Sends complex long double eigen vector data
//
// @param x			Complex eigen vector to send
// @param count	Number of vector elements
// @param id		Processor ID
// @param MTAG		Message tag
//**************************************************************************
void mpiProcess::send(const VectorCLD& x, int count, int id, int MTAG){
#ifdef HAVE_MPI
	MPI_Send(x.data(), count, MPI::LONG_DOUBLE_COMPLEX, id, MTAG, MPI_COMM_WORLD);
#endif
}

//**************************************************************************
// Receives complex long double eigen vector data
//
// @param x			Complex eigen vector to send
// @param count	Number of vector elements
// @param id		Processor ID
// @param MTAG		Message tag
//
// MPI requires me to pass in the actual vector to recieve the data. Thats
// why the data that you recieve is an input argument to the method. If 
// HAVE_MPI is not defined then the method just returns an empty vector.
//**************************************************************************
VectorCLD mpiProcess::recv(const VectorCLD& x, int count, int islave, int MTAG){
	// Inits a zero return vector
	VectorCLD xRet;
#ifdef HAVE_MPI
	MPI_Status status;
	MPI_Recv(x.data(), count, MPI::LONG_DOUBLE_COMPLEX, islave, MTAG, 
		MPI_COMM_WORLD, &status);
	xRet = x;
#endif
	return xRet;
}
