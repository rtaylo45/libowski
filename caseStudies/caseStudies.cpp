#include "mpiProcess.h"
#include "modelMesh.h"
#include "speciesDriver.h"
#include "meshCellData.h"
#include "species.h"

//**************************************************************************
// A "Lump" depletion case with a single cell. 
//
//	Problem equations:
//		dCi/dt = sum^{J}_{j=1} A_{i,j}*C_{j}
//
//	Domaine:
//			x = [0, 220]	cm
//			y = [0, 220]	cm
//			t = [0, 10]		y
//
//	Initial conditions and BC's:
//			C_{H-1}		= 2.74768E28
//			C_{B-10}		= 1.001E25
//			C_{O-16}		= 2.7503E28
//			C_{Zr-91}	= 4.14467E26
//			C_{Zr-96}	= 1.03432E26
//			C_{U-235}	= 1.909E26
//			C_{U-238}	= 6.592E27
//
//	Solution:
//		Calculated by matlab using year long time steps
// 
//**************************************************************************
void singleCellDepletion(int myid, std::string solverType){
	double t = 0.0;
	int steps = 10;
	double depletionTime = 10.*365.; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps;
	int xCells = 1, yCells = 1;
	double xLength = 220.0/1000., yLength = 220.0/1000.;
	std::vector<int> ids;
	std::string path = getDataPath();
	std::string outputFileName = "caseStudy1.out";
	std::string outputFileNameMatlab = "caseStudy1"+solverType+"Substeps12.csv";
	std::ofstream outputFile, outputFileMatlab;
	outputFile.open(outputFileName, std::ios_base::app);
	outputFileMatlab.open(outputFileNameMatlab);
	const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

	// Sets file paths for the input data	
	std::string speciesNamesFile = path + "speciesInputNames.dat";
	std::string speciesDecayFile = path + "speciesInputDecay.dat";
	std::string speciesTransFile = path + "speciesInputTrans.dat";

	// Builds the model
	modelMesh model(xCells, yCells, xLength, yLength);

	// Builds the species object
	speciesDriver spec = speciesDriver(&model);

	// Sets the neutron flux
	model.setSystemNeutronFlux(1.e13);

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);

	// Adds the speices
	ids = spec.addSpeciesFromFile(speciesNamesFile);
	// Sets the species sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);

	outputFile << "solverName: " << solverType << std::endl;
	outputFile.precision(16); 

	// Writes transition matrix and initial condition
	spec.writeTransitionMatrixToFile("transitionMatrixCaseStudy1.csv");
	spec.writeInitialConditionToFile("initialConditionCaseStudy1.csv");
	MatrixD solData = MatrixD::Zero(ids.size()+1,10);

	// Loops over the time steps to solve
	for (int k = 0; k < steps; k++){
		t = t + dt;
		spec.solve(t);
		if (myid == 0){
			std::cout << "solved" << std::endl;
			outputFile << "Time: " << t << std::endl;
			std::cout << "Time: " << t << std::endl;
			for (int id = 0; id < ids.size(); id++){
				std::string name = spec.getSpeciesName(0, 0, ids[id]);
				double con = spec.getSpecies(0, 0, ids[id]);
				outputFile << name << " " << con << std::endl;
				solData(id, k) = con;
			}
		}
	}
	outputFile << " " << std::endl;
	outputFileMatlab << solData.format(CSVFormat);
}

//**************************************************************************
// A pip depletion case
//
//	Problem equations:
//		dCi/dt = -v*dCi/dx + sum^{J}_{j=1} A_{i,j}*C_{j}
//
//	Domaine:
//			x = [0, 4]		m
//			t = [0, 10]		y
//
//	Initial conditions and BC's:
//			C_{H-1}		= 2.74768E28
//			C_{B-10}		= 1.001E25
//			C_{O-16}		= 2.7503E28
//			C_{Zr-91}	= 4.14467E26
//			C_{Zr-96}	= 1.03432E26
//			C_{U-235}	= 1.909E26
//			C_{U-238}	= 6.592E27
//
//	Solution:
//		Calculated by matlab using year long time steps
// 
//**************************************************************************
void pipeDepletion(int myid, std::string solverType){
	double t = 0.0;
	int steps = 1;
	double depletionTime = 100.; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps;
	int xCells = 1, yCells = 10;
	double velocity = 0.5;
	double xLength = 1.0, yLength = 6.0;
	std::vector<int> ids;

	// Files for the source terms and species names	
	std::string speciesNamesFile = getDataPath() + "speciesInputNames.dat";
	std::string speciesDecayFile = getDataPath() + "speciesInputDecay.dat";
	std::string speciesTransFile = getDataPath() + "speciesInputTrans.dat";

	// Builds the model mesh
	modelMesh model(xCells, yCells, xLength, yLength);

	// Inits the species driver
	speciesDriver spec = speciesDriver(&model);

	// Species IDs
	ids = spec.addSpeciesFromFile(speciesNamesFile);

	// Adds periodic boundary condions to all isotopes
	spec.setBoundaryCondition("periodic","south", ids);
	spec.setBoundaryCondition("periodic","north", ids);

	// Sets all the decay and transmutation sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);

	// Sets the neutron flux
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			meshCell* cell = model.getCellByLoc(i,j);
			double y = cell->y;
			double y1 = y - model.dy/2.;
			double y2 = y + model.dy/2.;
			double s = (1./model.dy)*(yLength/M_PI)*(cos(M_PI*y1/yLength) - 
			cos(M_PI*y2/yLength));
			model.setCellNeutronFlux(i, j, 1.e13*s);
		}
	}

	// set y velocity
	model.setConstantYVelocity(velocity);

	// Loops to solve the problem
	for (int k = 0; k < steps; k++){
		t = t + dt;
		if (myid == 0){std::cout << t << std::endl;};
		spec.solve(t);
	}
	if (myid == 0){
		// Loops to print results
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				std::cout << i << " " << j << std::endl;
				// Loops over the species 
				for (int id = 0; id < ids.size(); id++){
					std::string name = spec.getSpeciesName(i, j, ids[id]);
					double con = spec.getSpecies(i, j, ids[id]);
					std::cout << name << " " << con << std::endl;
				}
			}
		}
	}
}

//*****************************************************************************
// Main test
//*****************************************************************************
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;
	//std::vector<std::string> solvers {"taylor"};
	//std::vector<std::string> solvers {"CRAM", "hyperbolic", "parabolic", "pade-method1",
	//"pade-method2"};
	std::vector<std::string> solvers {"CRAM"};

	// Loops over different solvers
	for (std::string &solverType : solvers){
		if(myid == 0){std::cout << solverType << std::endl;};
		//singleCellDepletion(myid, solverType);
		pipeDepletion(myid, solverType);
	}

	mpi.finalize();
}
