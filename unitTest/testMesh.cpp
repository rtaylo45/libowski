#include <assert.h>
#include <iostream>
#include <string>
#include <vector>

#include "mpiProcess.h"
#include "modelMesh.h"
#include "speciesDriver.h"
#include "meshCell.h"
#include "species.h"
#include "sys.h"
#include "matrixTypes.h"

using namespace Eigen;

//*****************************************************************************
// test the init of the mesh class. For now i kinda just test to make sure
// the function run.
//*****************************************************************************
void testInit(){

	modelMesh mesh(2, 5, 1.0, 1.0);
	mesh.setConstantXVelocity(2.0, 0);
	mesh.setConstantXVelocity(2.0, 1);
	mesh.setConstantXVelocity(2.0);
	mesh.setConstantYVelocity(0.1);
	mesh.setConstantYVelocity(0.2, 1);
	mesh.clean();

}

//*****************************************************************************
// test the init of the species driver. Test to make sure multiple species
// can be added. initial concentrations are added correctly. Molar masses
// are added currect and sources are set right. 
//*****************************************************************************
void testSpeciesDriver(){
	int xCells = 2, yCells = 5;
	double xLength = 1.0, yLength = 1.0;
	double spec1InitCon = 1.0, spec2InitCon = 2.0;
	double spec1MM = 2.0, spec2MM = 3.0;
	int specID1, specID2;
	double spec1Con, spec2Con;
	std::vector<double> spec1Coeffs = {5.0, 10.0};
	std::vector<double> spec2Coeffs = {6.0, 11.0};
	double spec1S = 50.0, spec2S = 100.0;

	modelMesh model(xCells, yCells, xLength, yLength);
	speciesDriver spec = speciesDriver(&model);
	specID1 = spec.addSpecies(spec1MM, spec1InitCon, 0.0, "spec1");
	specID2 = spec.addSpecies(spec2MM, spec2InitCon, 0.0, "spec2");

	// Sets sourse terms for the model
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			spec.setSpeciesSource(i,j,specID1,	spec1Coeffs, spec1S);
			spec.setSpeciesSource(i,j,specID2,	spec2Coeffs, spec2S);
		}
	}
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			// Gets species pointers
			species* spec1 = spec.getSpeciesPtr(i, j, specID1);
			species* spec2 = spec.getSpeciesPtr(i, j, specID2);

			// Gets species Concentrations
			spec1Con = spec.getSpecies(i, j, specID1);;
			spec2Con = spec.getSpecies(i, j, specID2);;

			// Makes sure all species concentrations are right
			assert(1.0 == spec1Con);
			assert(2.0 == spec2Con);
			// Makes sure all species molar masses are right
			assert(spec1->MM == spec1MM);
			assert(spec2->MM == spec2MM);
		}
	}
	assert(specID1 == spec.getSpeciesID("spec1"));
	assert(specID2 == spec.getSpeciesID("spec2"));
	
	model.clean();
	spec.clean();

}

//*****************************************************************************
// Test adding species from file
//*****************************************************************************
void testAddSpeciesFromFile(){
	int xCells = 1, yCells = 1;
	double xLength = 1.0, yLength = 1.0;
	std::string data = getDataPath() + "test/";
	std::vector<int> ids; std::vector<std::vector<double>> vectDecay;
	std::vector<std::vector<double>> vectTran;
	std::vector<double> decayVector;
	std::vector<double> tranVector;
	std::vector<std::string> testName = {"U-235", "Xe-135m", "Xe-135", "I-135"};
	std::vector<double> testMM = {235.0439, 134.9072, 134.9072, 134.91};
	// Decay vectors
	std::vector<double> U235Decay = {-3.1209e-17, 0.0, 0.0, 0.0};
	std::vector<double> Xe135mDecay = {0.0, -0.00075553, 0.0, 4.8381e-06};
	std::vector<double> Xe135Decay = {0.0, 0.00075326, -2.1066e-05, 2.4468e-05};
	std::vector<double> I135Decay = {0.0, 0.0, 0.0, -2.9306e-05};
	vectDecay = {U235Decay, Xe135mDecay, Xe135Decay, I135Decay};
	// Transition vectors
	std::vector<double> U235Tran = {0.0, 0.0, 0.0, 0.0};
	std::vector<double> Xe135mTran = {4.9184e-26, 0.0, 7.3058e-26, 0.0};
	std::vector<double> Xe135Tran = {9.6002e-26, 0.0, 0.0, 0.0};
	std::vector<double> I135Tran = {8.8728e-25, 0.0, 1.9667e-31, 0.0};
	vectTran = {U235Tran, Xe135mTran, Xe135Tran, I135Tran};
	
	
	std::string speciesNamesFile = data + "speciesInputNamesSmall.txt";
	std::string speciesDecayFile = data + "speciesInputDecaySmall.txt";
	std::string speciesTransFile = data + "speciesInputTransSmall.txt";

	modelMesh model(xCells, yCells, xLength, yLength);
	speciesDriver spec = speciesDriver(&model);

	model.setSystemNeutronFlux(1.0);

	ids = spec.addSpeciesFromFile(speciesNamesFile);
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);

	for (int i = 0; i < ids.size(); i++){
		species* thisSpec = spec.getSpeciesPtr(0, 0, ids[i]);
		meshCell* cell = model.getCellByLoc(0, 0);
		assert(thisSpec->name == testName[i]);
		assert(thisSpec->MM == testMM[i]);
		decayVector = vectDecay[i];
		tranVector = vectTran[i];
		for (int otherID = 0; otherID < ids.size(); otherID++){
			double testDecayCoeff = thisSpec->getTransitionCoeff(otherID, 0, 
				cell->getScalarData());
			double testTranCoeff = thisSpec->getTransitionCoeff(otherID, 1, 
				cell->getScalarData());
			assert(testDecayCoeff == decayVector[otherID]);
			assert(testTranCoeff == tranVector[otherID]);
		}
	}
}

//*****************************************************************************
// Test setting up gas source term from file
//*****************************************************************************
void testSetGasSpargingAndWallDepositionFromFile(){
	int xCells = 1, yCells = 1;
	double xLength = 1.0, yLength = 1.0;
	double intarea = 1./10., temp = 100., voidFract = 0.1, surfarea = 10.;
	double k1 = 2.0, k2 = 1.0, h1 = 5.0e-1, h2 = 6.0e-2;
	double liqCoeff1 = k1*intarea/(1.-voidFract);
	double liqCoeff2 = k2*intarea/(1.-voidFract);
	double gasCoeff1 = k1*intarea*h1*idealGasR*temp/voidFract;
	double gasCoeff2 = k2*intarea*h2*idealGasR*temp/voidFract;
	std::vector<int> ids;
	std::string data = getDataPath() + "test/";
	std::vector<double> testGas = {-liqCoeff2, 0, 0, gasCoeff2, 0, -liqCoeff1, 
		gasCoeff1, 0, 0, liqCoeff1, -gasCoeff1, 0, liqCoeff2, 0, 0, -gasCoeff2};
	std::vector<double> testWall = {-20, 0, 20, 0, 0, -50, 0, 50, 20, 0, -20,
		0, 0, 50, 0, -50};
	
	std::string speciesNamesFile = data + "speciesInputNamesSmall.txt";

	modelMesh model(xCells, yCells, xLength, yLength);
	speciesDriver spec = speciesDriver(&model);

	model.setSystemGasInterfacialAreaCon(intarea);
	model.setSystemTemperature(temp);
	model.setSystemGasVoidFraction(voidFract);
	model.setSystemWallInterfacialAreaCon(surfarea);

	ids = spec.addSpeciesFromFile(speciesNamesFile);
	spec.setGasSpargingFromFile(data + "speciesInputGasSpargingSmall.txt");
	spec.setWallDepositionFromFile(data + "speciesInputWallDepositionSmall.txt");

	int index = 0;
	for (int i = 0; i < ids.size(); i++){
		species* thisSpec = spec.getSpeciesPtr(0, 0, ids[i]);
		meshCell* cell = model.getCellByLoc(0, 0);
		for (int otherID = 0; otherID < ids.size(); otherID++){
			double gasCoeff = thisSpec->getTransitionCoeff(otherID, 0, 
				cell->getScalarData());
			double wallCoeff = thisSpec->getTransitionCoeff(otherID, 1, 
				cell->getScalarData());
			if (gasCoeff != 0){assert(isApprox(gasCoeff, testGas[index]));};
			if (wallCoeff != 0){assert(isApprox(wallCoeff, testWall[index]));};
			index ++;
		}
	}
}

//*****************************************************************************
// Main test
//*****************************************************************************
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;

	testInit();
	testSpeciesDriver();
	testAddSpeciesFromFile();
	testSetGasSpargingAndWallDepositionFromFile();

	mpi.finalize();
}
