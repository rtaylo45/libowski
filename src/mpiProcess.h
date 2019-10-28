//*****************************************************************************
// Author: Zack Taylor
//
// MPI wrapper class: Used to make using MPI easier and for easier ability 
// to run code in serial or parallel. Only code that is inside #ifdef and 
// #endif will be scene by the compiler and exicuted if the HAVE_MPI compiler
// variable is set in the cmake build process. HAVE_MPI needs to be set in 
// the CMakeList.txt file as a compile definition. 
//
// At the end of the file we define the mpi object with 
//
// mpiProcess mpi;
//
// This is the object that needs to be used by all other things in the code.
//*****************************************************************************
#ifndef MPIPROCESS_H
#define MPIPROCESS_H
#include <complex>
#include <Eigen/Sparse>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace Eigen;

class mpiProcess{

	// Public attributes
	public:
	// logical for initilization
	bool isInit = false;
	// Number of total processors
	int size = 1;
	// ID of the current processor
	int rank = 0;

	//**************************************************************************
	// Constructor
	//**************************************************************************
	public:
	mpiProcess();

	// Definition of methods
	public:
	//**************************************************************************
	// Finalize MPI
	//**************************************************************************
	void finalize();
	//**************************************************************************
	// Sends complex eigen vector data
	//**************************************************************************
	void send(VectorXcd, int, int, int);
	//**************************************************************************
	// Receives complex eigen vector data
	//**************************************************************************
	VectorXcd recv(VectorXcd, int, int, int);

	//**************************************************************************
	// Sends complex eigen sparse matrix data
	//**************************************************************************
	void send(SparseMatrix<std::complex<double>>, int, int, int);
	//**************************************************************************
	// Receives complex eigen sparse matrix data
	//**************************************************************************
	SparseMatrix<std::complex<double>> recv(SparseMatrix<std::complex<double>>, 
		int, int, int);

	//**************************************************************************
	// Initialization of the mpi object
	//**************************************************************************
	private:
	void initMPI();
};
// points to the global mpi variable in mpiProcess.cpp
extern mpiProcess mpi;
#endif
