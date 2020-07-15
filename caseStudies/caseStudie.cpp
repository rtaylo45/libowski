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

//*****************************************************************************
// Main test
//*****************************************************************************
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	singleCellDepletion(myid);

	mpi.finalize();
}
