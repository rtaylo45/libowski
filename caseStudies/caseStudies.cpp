#include "mpiProcess.h"
#include "modelMesh.h"
#include "speciesDriver.h"
#include "meshCellData.h"
#include "species.h"

void singleCellDepletion(int myid){
	double t = 0.0;
	int steps = 10;
	double depletionTime = 100.; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps;
	int xCells = 1, yCells = 1;
	double xLength = 1.0, yLength = 1.0;
	std::vector<int> ids;
	
	//std::string speciesNamesFile = getDataPath() + "speciesInputNamesSmall.txt";
	//std::string speciesDecayFile = getDataPath() + "speciesInputDecaySmall.txt";
	//std::string speciesTransFile = getDataPath() + "speciesInputTransSmall.txt";
	std::string speciesNamesFile = getDataPath() + "speciesInputNames.dat";
	std::string speciesDecayFile = getDataPath() + "speciesInputDecay.dat";
	std::string speciesTransFile = getDataPath() + "speciesInputTrans.dat";

	modelMesh model(xCells, yCells, xLength, yLength);
	speciesDriver spec = speciesDriver(&model);

	ids = spec.addSpeciesFromFile(speciesNamesFile);
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);

	// Sets the neutron flux
	model.setSystemNeutronFlux(1.e13);

	for (int k = 0; k < steps; k++){
		t = t + dt;
		spec.solve(t);
	}
	for (int id = 0; id < ids.size(); id++){
		std::string name = spec.getSpeciesName(0, 0, ids[id]);
		double con = spec.getSpecies(0, 0, ids[id]);
		std::cout << name << " " << con << std::endl;
	}
}

void pipeDepletion(int myid){
	double t = 0.0;
	int steps = 5;
	double depletionTime = 0.05; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps;
	int xCells = 1, yCells = 50;
	double velocity = 0.5;
	double xLength = 1.0, yLength = 100.0;
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

	//singleCellDepletion(myid);
	pipeDepletion(myid);

	mpi.finalize();
}
