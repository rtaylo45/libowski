#include "mpiProcess.h"

// Global mpi object
mpiProcess mpi;

//*****************************************************************************
// Constructor
//*****************************************************************************
mpiProcess::mpiProcess(){initMPI();};

//*****************************************************************************
// Init the mpi process
//*****************************************************************************
void mpiProcess::initMPI(){
#ifdef HAVE_MPI
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	isInit = true;
#endif
}

//*****************************************************************************
// Finalize the mpi process
//*****************************************************************************
void mpiProcess::finalize(){
#ifdef HAVE_MPI
	MPI_Finalize();
	isInit = false;
#endif
}

//*****************************************************************************
// Sends complex long double eigen sparse matrix data
//
// @param matIn	Complex eigen sparse matrix to send
// @param id		Processor ID
// @param tagBase	Base for the mpi tage
//*****************************************************************************
void mpiProcess::send(SparseMatrixCLD& matIn, int id, int tagBase){
#ifdef HAVE_MPI
	MPI_Request req;
	MPI_Status status;
	matIn.makeCompressed();
	int mtag0 = tagBase, mtag1 = tagBase+1, mtag2 = tagBase+2, mtag3 = tagBase+3;
	int rows = matIn.rows(), cols = matIn.cols(), nnz = matIn.nonZeros();
	assert(rows==matIn.innerSize() && cols==matIn.outerSize());
	assert(matIn.outerIndexPtr()[cols]==nnz);
	int shape[3] = {rows, cols, nnz};

	MPI_Isend(matIn.innerIndexPtr(),nnz, MPI::INT, id, mtag2, MPI_COMM_WORLD, &req);
	MPI_Isend(shape,					 3, MPI::INT, id, mtag0, MPI_COMM_WORLD, &req);
	MPI_Isend(matIn.outerIndexPtr(),cols, MPI::INT, id, mtag3, MPI_COMM_WORLD, &req);
	MPI_Isend(matIn.valuePtr(),		 nnz, MPI::LONG_DOUBLE_COMPLEX, id, mtag1, MPI_COMM_WORLD, &req);
	MPI_Wait(&req, &status);
#endif
}

//*****************************************************************************
// Receives complex long double eigen sparse matrix data
//
// @param matOut	Complex eigen sparse matrix to receive
// @param id		Processor ID
// @param tagBase	Base for the mpi tage
//
//*****************************************************************************
void mpiProcess::recv(SparseMatrixCLD& matOut, int islave, int tagBase){
#ifdef HAVE_MPI
	MPI_Status status;
	int shape[3];
	MPI_Recv(shape, 3, MPI_INT, islave, tagBase, MPI_COMM_WORLD, &status);
	int rows = shape[0], cols = shape[1], nnz = shape[2];
	matOut.resize(rows, cols);
	matOut.reserve(nnz);
	MPI_Recv(matOut.valuePtr(), nnz, MPI::LONG_DOUBLE_COMPLEX, islave, tagBase+1, 
		MPI_COMM_WORLD, &status);
	MPI_Recv(matOut.innerIndexPtr(), nnz, MPI::INT, islave, tagBase+2, 
		MPI_COMM_WORLD, &status);
	MPI_Recv(matOut.outerIndexPtr(), cols, MPI::INT, islave, tagBase+3, 
		MPI_COMM_WORLD, &status);
	matOut.outerIndexPtr()[cols] = nnz;
#endif
}

//*****************************************************************************
// Sends complex long double eigen vector data
//
// @param vectIn	Complex eigen vector to send
// @param id		Processor ID
// @param tagBase	Base for the mpi tage
//
//*****************************************************************************
void mpiProcess::send(VectorCLD& vectIn, int id, int tagBase){
#ifdef HAVE_MPI
	int rows = vectIn.rows();
	int shape[1] = {rows};
	MPI_Send(shape, 1, MPI::INT, id, tagBase, MPI_COMM_WORLD);
	MPI_Send(vectIn.data(), rows, MPI::LONG_DOUBLE_COMPLEX, id, tagBase+1, 
		MPI_COMM_WORLD);
#endif
}

//*****************************************************************************
// Receives complex long double eigen vector data
//
// @param out		Complex eigen vector to receive
// @param id		Processor ID
// @param tagBase	Base for the mpi tage
//
//*****************************************************************************
void mpiProcess::recv(VectorCLD& out, int islave, int tagBase){
#ifdef HAVE_MPI
	MPI_Status status;
	int shape[1];
	MPI_Recv(shape, 1, MPI::INT, islave, tagBase, MPI_COMM_WORLD, &status);
	int rows = shape[0];
	out.resize(rows);
	MPI_Recv(out.data(), rows, MPI::LONG_DOUBLE_COMPLEX, islave, tagBase+1, 
		MPI_COMM_WORLD, &status);
#endif
}
