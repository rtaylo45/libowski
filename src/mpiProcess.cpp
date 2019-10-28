#include <iostream>
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
// Sends complex eigen sparse matrix data
//
// @param x			Complex eigen sparse matrix to send
// @param count	Number of matrix elements
// @param id		Processor ID
// @param MTAG		Message tag
//**************************************************************************
void mpiProcess::send(SparseMatrix<std::complex<double>> x, int count, 
	int id, int MTAG){
#ifdef HAVE_MPI
	MPI_Send(x.data(), count, MPI::DOUBLE_COMPLEX, id, MTAG, MPI_COMM_WORLD);
#endif
}

//**************************************************************************
// Receives complex eigen sparse matrix data
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
SparseMatrix<std::complex<double>> mpiProcess::recv(SparseMatrix<std::complex<double>> 
	x, int count, int islave, int MTAG){
	// Inits a zero return vector
	SparseMatrix<std::complex<double>> xRet;
#ifdef HAVE_MPI
	MPI_Status status;
	MPI_Recv(x.data(), count, MPI::DOUBLE_COMPLEX, islave, MTAG, 
		MPI_COMM_WORLD, &status);
	xRet = x;
#endif
	return xRet;
}

//**************************************************************************
// Sends complex eigen vector data
//
// @param x			Complex eigen vector to send
// @param count	Number of vector elements
// @param id		Processor ID
// @param MTAG		Message tag
//**************************************************************************
void mpiProcess::send(VectorXcd x, int count, int id, int MTAG){
#ifdef HAVE_MPI
	MPI_Send(x.data(), count, MPI::DOUBLE_COMPLEX, id, MTAG, MPI_COMM_WORLD);
#endif
}

//**************************************************************************
// Receives complex eigen vector data
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
VectorXcd mpiProcess::recv(VectorXcd x, int count, int islave, int MTAG){
	// Inits a zero return vector
	VectorXcd xRet;
#ifdef HAVE_MPI
	MPI_Status status;
	MPI_Recv(x.data(), count, MPI::DOUBLE_COMPLEX, islave, MTAG, 
		MPI_COMM_WORLD, &status);
	xRet = x;
#endif
	return xRet;
}
