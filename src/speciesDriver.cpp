//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include "speciesDriver.h"

//*****************************************************************************
// Constructor
//
// @param model A pointer to the mesh model
//*****************************************************************************
speciesDriver::speciesDriver(modelMesh* model){
	modelPtr = model;
	expSolver = matrixExponentialFactory::getExpSolver("CRAM");
	intSolver = integratorFactory::getIntegrator("implicit", "BDF2");
}

//*****************************************************************************
// Sets the matrix exponential solver
//
// @param solverName	The name of the materix exponential solver type
// @param krylovFlag	Bool used to set if the krylov subspace is to be used
// @param krylovDim	Dimension of the subspace
//*****************************************************************************
void speciesDriver::setMatrixExpSolver(std::string solverName, bool krylovFlag,
	int krylovDim){
	expSolver = matrixExponentialFactory::getExpSolver(solverName, krylovFlag,
		krylovDim);
}

//*****************************************************************************
// Sets the krylov subsapce dimension of the matrix exp solver
//
// @param dim	Dimension of the krylov subspace
//*****************************************************************************
void speciesDriver::setKrylovSubspaceDimension(int dim){
	expSolver->setKrylovSubspaceDimension(dim);
}

//*****************************************************************************
// Sets the integrator solver
//
// @param method		The method of the solver
// @param solverName	The name of the materix exponential solver type
//*****************************************************************************
void speciesDriver::setIntegratorSolver(std::string method, std::string 
	solverName){
	intSolver = integratorFactory::getIntegrator(method, solverName);
}

//*****************************************************************************
// Writes out the base line transition matrix to a csv file. This matrix
// is not multiplied by the time step size. This function must be called
// After all of the speices source terms are set and all mesh parameters 
// are set. Call this right before you call the solve method
//
// @param fname	The name of the file to write the matrix to. This file will
//						be in the csv format
//*****************************************************************************
void speciesDriver::writeTransitionMatrixToFile(const std::string fname){
	MatrixD dA;
	bool augmented = true;

	if (mpi.rank == 0){
		A = buildTransMatrix(augmented, 0.0);
		dA = Eigen::MatrixXd(A);
		writeCSV(dA, fname);
	}	
}
//*****************************************************************************
// Writes out the Initial condition vector
//
// @param fname	The name of the file to write the matrix to. This file will
//						be in the csv format
//*****************************************************************************
void speciesDriver::writeInitialConditionToFile(const std::string fname){
	VectorD v;
	bool augmented = true;

	if (mpi.rank == 0){
		v = buildInitialConditionVector(augmented);
		writeCSV(v, fname);
	}	
}

//*****************************************************************************
// Sets the flux limiter function
//
// @param limiterName	The name of the flux limiter
//*****************************************************************************
void speciesDriver::setFluxLimiter(std::string limiterName){
	int limiterID = -1;	
	if (limiterName == "superbee"){ limiterID = 0;};
	if (limiterName == "VanLeer"){ limiterID = 1;};
	if (limiterName == "Van Albada"){ limiterID = 2;};
	if (limiterName == "Min-Mod"){ limiterID = 3;};
	if (limiterName == "Sweby"){ limiterID = 4;};
	if (limiterName == "First order upwind"){ limiterID = 5;};
	if (limiterName == "muscl"){ limiterID = 6;};

	fluxLim.setLimiterFunction(limiterID);
}

//*****************************************************************************
// Adds a species to the model
//
// @param molarMass  Molar mass of species [g/mol]
// @param [initCon]  Initial concentration [kg/m^3]
// @param [diffCoef]	Diffusion coefficient [m^2/s]
// @param name			Name of the species
// @param Transport	bool to set if the speices is to be transport with the 
//							velocity field
//*****************************************************************************
int speciesDriver::addSpecies(double molarMass, double initCon,
	double diffCoeff, std::string name, bool transport){
	// the Species ID
   int specID = numOfSpecs;

	// Generates default name if none is given
	if (name == "None"){
		name = "spec" + std::to_string(specID);
	}
   for (int i = 0; i < modelPtr->numOfxCells; i++){
      for (int j = 0; j < modelPtr->numOfyCells; j++){
         meshCell* cell = modelPtr->getCellByLoc(i,j);

         cell->addSpecies(molarMass, initCon, diffCoeff, name, transport);
      }
   }

	// Map the name to the spec ID
	specNameToID.insert(std::pair<std::string, int>(name, specID));
	// Mape the spec ID to the name
	specIDToName.insert(std::pair<int, std::string>(specID, name));
	// Makes sure that the species name does not already exist. This would make
	// the count 2.
	assert(specNameToID.count(name) == 1);
	assert(specIDToName.count(specID) == 1);

   numOfSpecs++;
   return specID;
}

//*****************************************************************************
// Add species from a file generated from pyLibowski. Returns a vector of species
// IDs
//
// @param fname	File location
//*****************************************************************************
std::vector<int> speciesDriver::addSpeciesFromFile(std::string fname){
   std::ifstream infile(fname);
	std::vector<int> specIDs;
   double mm, initCon, D;
   std::string name;

	// Checks to see if the file exists. throws error if not
	checkFileExists(fname);
		
	// Loop though lines
   for (std::string line; getline(infile, line);){
      std::istringstream iss(line);
      std::vector<std::string> result;
		// Split line and loop though it
      for (std::string s; iss >> s;){
         result.push_back(s);
      }
		// Skips comments
		if (result.at(0) != "#"){
			name = result.at(0); mm = stod(result.at(1)); initCon = stod(result.at(2)); 
			D = stod(result.at(3));
			// Adds the species
			specIDs.push_back(addSpecies(mm, initCon, D, name, true));
		}
   }
	return specIDs;
}

//*****************************************************************************
// Returns a pointer to the spcies object in a cell
//
// @param i       x index
// @param j       y index
// @param specID  Species ID
//*****************************************************************************
species* speciesDriver::getSpeciesPtr(int i, int j, int specID){
   meshCell* cell = modelPtr->getCellByLoc(i,j);
   species* specPtr = cell->getSpecies(specID);
   return specPtr;
}

//*****************************************************************************
// Gets the species concentration
//
// @param i       x index
// @param j       y index
// @param specID  Species ID
//*****************************************************************************
double speciesDriver::getSpecies(int i, int j, int specID){
   species* spec = getSpeciesPtr(i, j, specID);
   double specCon = spec->c;
   return specCon;
}

//*****************************************************************************
// Gets the species ID for a given species name. If there is not a match
// then it throws an error
//
// @param name		Species name
//*****************************************************************************
int speciesDriver::getSpeciesID(std::string name){
	// Makes sure the name is in the map
	assert(specNameToID.count(name) == 1);
	int id = specNameToID[name];
	return id;
}

//*****************************************************************************
// Sets the species concentration
//
// @param i       x index
// @param j       y index
// @param specID  Species ID
// @param specCon	Concentration [kg/m^3]
//*****************************************************************************
void speciesDriver::setSpeciesCon(int i, int j, int specID, double specCon){
   species* spec = getSpeciesPtr(i, j, specID);
   spec->c = specCon;
}

//*****************************************************************************
// Gets the species name
//
// @param specID  Species ID
//*****************************************************************************
std::string speciesDriver::getSpeciesName(int specID){
	// Makes sure the name is in the map
	assert(specIDToName.count(specID) == 1);
	std::string name = specIDToName[specID];
   return name;
}

//*****************************************************************************
// Sets general source terms for a species in a cell
//
// @param i					x index
// @param j					y index
// @param specID			Species ID
// @param coeffs			An array of species source coefficients
//								[1/s]
// @param s					Constant source in cell [kg/m^3/s]
//*****************************************************************************
void speciesDriver::setSpeciesSource(int i, int j, int specID, std::vector<double>
      coeffs, double s){
   assert(coeffs.size() == numOfSpecs);
   species* spec = getSpeciesPtr(i, j, specID);
	spec->addGenericSourceTerm(coeffs);
   spec->s = s;
}

//*****************************************************************************
// Internal function that sets the decay source terms for a species in a cell
//
// The name is passed in so that an assertion is made that the species name 
// equals the name in the file
//
// @param i					x index
// @param j					y index
// @param specID			Species ID
// @param name				Species name
// @param coeffs			An array of species source coefficients
//								[1/s]
//*****************************************************************************
void speciesDriver::setDecaySource(int i, int j, int specID, std::string name,
		std::vector<double> coeffs){
   species* spec = getSpeciesPtr(i, j, specID);
	assert(spec->name == name);
	assert(numOfSpecs == coeffs.size());
   spec->addGenericSourceTerm(coeffs);
}

//*****************************************************************************
// Internal function that sets the trans source terms for a species in a cell
//
// The name is passed in so that an assertion is made that the species name 
// equals the name in the file
//
// @param i					x index
// @param j					y index
// @param specID			Species ID
// @param name				Species name
// @param coeffs			A vector of species source coefficients
//								[cm^2]
//*****************************************************************************
void speciesDriver::setTransSource(int i, int j, int specID, std::string name,
		std::vector<double> coeffs){
   species* spec = getSpeciesPtr(i, j, specID);
	assert(spec->name == name);
	assert(numOfSpecs == coeffs.size());
   spec->addNIRSourceTerm(coeffs);
}

//*****************************************************************************
// Function that sets the wall deposition model for a single cell
//
// @param i					x index
// @param j					y index
// @param coeffs			A vector of mass transfer coefficients [m/s]
// @param liquidIDs		Vector of liquid IDs
// @param surfaceIDs		Vector of surface IDs must be the same order as 
//								liquid IDs
//	@param infSink			Bools for infinite sink assumption [false]
//*****************************************************************************
void speciesDriver::setWallDeposition(int i, int j, std::vector<double> coeffs,
		std::vector<int> liquidIDs, std::vector<int> surfaceIDs, 
		std::vector<bool> infSinks){
	std::vector<bool> infSinks_;
	assert(liquidIDs.size() == surfaceIDs.size());
	assert(liquidIDs.size() == coeffs.size());
	if (infSinks.size() != 0){
		assert(infSinks.size() == liquidIDs.size());
		infSinks_ = infSinks;
	}
	else {
		for (int i = 0; i < coeffs.size(); i++){infSinks_.push_back(false);};
	}

	for (int index = 0; index < coeffs.size(); index++){
		species* specLiq = getSpeciesPtr(i, j, liquidIDs[index]);
		species* specSurf = getSpeciesPtr(i, j, surfaceIDs[index]);
		specLiq->addWallDepositionSourceTerm(coeffs[index], liquidIDs[index], 
			liquidIDs[index], surfaceIDs[index], infSinks_[index]);
		specSurf->addWallDepositionSourceTerm(coeffs[index], surfaceIDs[index],
			liquidIDs[index], surfaceIDs[index], infSinks_[index]);
		specSurf->transport = false;
	}

}

//*****************************************************************************
// Function that sets the wall deposition model for the whole system
//
// @param coeffs			A vector of mass transfer coefficients [m/s]
// @param liquidIDs		Vector of liquid IDs
// @param surfaceIDs		Vector of surface IDs must be the same order as 
//								liquid IDs
//	@param infSinks		Bool for infinite sink assumption 
//*****************************************************************************
void speciesDriver::setWallDeposition(std::vector<double> coeffs,
		std::vector<int> liquidIDs, std::vector<int> surfaceIDs, std::vector<bool> 
		infSinks){

	// Loop over cells
	for (int i = 0; i < modelPtr->numOfxCells; i++){
		for (int j = 0; j < modelPtr->numOfyCells; j++){
			setWallDeposition(i, j, coeffs, liquidIDs, surfaceIDs, infSinks);
		}
	}

}

//*****************************************************************************
// Function that sets the gas sparging model for a single cell
//
// @param i					x index
// @param j					y index
// @param mCoeffs			A vector of mass transfer coefficients [m/s]
// @param HCoeffs			A vector of Henry laws coefficients		[mol/m^3/Pa]
// @param liquidIDs		Vector of liquid IDs
// @param gasIDs			Vector of gas IDs must be the same order as 
//								liquid IDs
//*****************************************************************************
void speciesDriver::setGasSparging(int i, int j, std::vector<double> mCoeffs,
		std::vector<double> HCoeffs, std::vector<int> liquidIDs, 
		std::vector<int> gasIDs){
	assert(liquidIDs.size() == gasIDs.size());
	assert(liquidIDs.size() == mCoeffs.size());
	assert(liquidIDs.size() == HCoeffs.size());

	for (int index = 0; index < mCoeffs.size(); index++){
		species* specLiq = getSpeciesPtr(i, j, liquidIDs[index]);
		species* specGas = getSpeciesPtr(i, j, gasIDs[index]);
		specLiq->addGasSpargingSourceTerm(mCoeffs[index], HCoeffs[index], liquidIDs[index], 
			liquidIDs[index], gasIDs[index]);
		specGas->addGasSpargingSourceTerm(mCoeffs[index], HCoeffs[index], gasIDs[index],
			liquidIDs[index], gasIDs[index]);
	}

}

//*****************************************************************************
// Function that sets the gas sparging model for the whole system
//
// @param mCoeffs			A vector of mass transfer coefficients [m/s]
// @param HCoeffs			A vector of Henry laws coefficients		[mol/m^3/Pa]
// @param liquidIDs		Vector of liquid IDs
// @param gasIDs			Vector of gas IDs must be the same order as 
//								liquid IDs
//*****************************************************************************
void speciesDriver::setGasSparging(std::vector<double> mCoeffs, 
		std::vector<double> HCoeffs, std::vector<int> liquidIDs, 
		std::vector<int> gasIDs){

	// Loop over cells
	for (int i = 0; i < modelPtr->numOfxCells; i++){
		for (int j = 0; j < modelPtr->numOfyCells; j++){
			setGasSparging(i, j, mCoeffs, HCoeffs, liquidIDs, gasIDs);
		}
	}

}

//*****************************************************************************
// Sets the source terms for a species in the entire problem domain
//
// @param decayfname		File location of the decay file
// @param transfname		File location of the trans file (optional)
//*****************************************************************************
void speciesDriver::setSpeciesSourceFromFile(std::string decayfname, std::string transfname){
	std::vector<std::string> fileVect;
	std::string name;
	int specID;
	// the trans file is not required and has the default name None
	if (transfname == "None"){ 
		// Checks to see if the file exists. throws error if not
		checkFileExists(decayfname);
		// Adds the file to the fname vector
		fileVect.push_back(decayfname);
	}
	// Both files are entered
	else{
		// Checks to see if the files exists. throws error if not
		checkFileExists(decayfname);
		checkFileExists(transfname);
		// Adds the file to the fname vector
		fileVect.push_back(decayfname);
		fileVect.push_back(transfname);
	}
	
	// Loops though the files	
	for(int findex = 0; findex < fileVect.size(); findex++){
		std::ifstream infile(fileVect[findex]);	
		// Loops thought the lines of the file
		for (std::string line; getline(infile, line);){
   	   std::istringstream iss(line);
   	   std::vector<double> v;
   	   int l = 0;
			// Loops though each line to split the string
   	   for (std::string s; iss >> s;){
				// specID is the first entry
   	      if (l == 0){specID = stoi(s);};
				// name is the second entry
   	      if (l == 1){name = s;};
				// the rest of the entries are the source terms
   	      if (l > 1){
   	         v.push_back(stod(s));
   	      }
   	      l++;
   	   }

			// Loop over cells
			for (int i = 0; i < modelPtr->numOfxCells; i++){
				for (int j = 0; j < modelPtr->numOfyCells; j++){
					if (findex == 0){ 
						setDecaySource(i, j, specID, name, v);
					}
					else{
						setTransSource(i, j, specID, name, v);
					}
				}
			}
   	}
	}
}
//*****************************************************************************
// Sets up gas sparging for species in the entire problem domain
//
// @param gasfname		File location of the gas transport file
//*****************************************************************************
void speciesDriver::setGasSpargingFromFile(std::string gasfname){
	int liqID, gasID;
	double k, H;
	std::vector<int> liqIDs, gasIDs;
	std::vector<double> mCoeffs, hCoeffs;
	// Checks to see if the files exists. throws error if not
	checkFileExists(gasfname);
	
	// Loops thought the lines of the file
	std::ifstream infile(gasfname);	
	for (std::string line; getline(infile, line);){
      std::istringstream iss(line);
      std::vector<double> v;
      int l = 0;
		// Loops though each line to split the string
      for (std::string s; iss >> s;){
			// liquid spec name is the first entry
         if (l == 0){liqID = getSpeciesID(s);};
			// gas spec name is the second entry
         if (l == 1){gasID = getSpeciesID(s);};
			// the third entry is the mass transfer coefficient
         if (l == 2){k = stod(s);};
			// the last entry is the hentrys law coefficient
			if (l == 3){H = stod(s);};
         l++;
      }
		liqIDs.push_back(liqID); gasIDs.push_back(gasID);
		mCoeffs.push_back(k); hCoeffs.push_back(H);

   }
	setGasSparging(mCoeffs, hCoeffs, liqIDs, gasIDs);
}

//*****************************************************************************
// Sets up wall deposition for species in the entire problem domain
//
// @param wallfname		File location of the gas transport file
//*****************************************************************************
void speciesDriver::setWallDepositionFromFile(std::string wallfname){
	int liqID, surfID, integer;
	double k, H;
	bool infsink;
	std::vector<int> liqIDs, surfIDs;
	std::vector<double> mCoeffs, hCoeffs;
	std::vector<bool> infSinks;
	// Checks to see if the files exists. throws error if not
	checkFileExists(wallfname);
	
	// Loops thought the lines of the file
	std::ifstream infile(wallfname);	
	for (std::string line; getline(infile, line);){
      std::istringstream iss(line);
      std::vector<double> v;
      int l = 0;
		// Loops though each line to split the string
      for (std::string s; iss >> s;){
			// liquid spec name is the first entry
         if (l == 0){liqID = getSpeciesID(s);};
			// surface spec name is the second entry
         if (l == 1){surfID = getSpeciesID(s);};
			// the third entry is the mass transfer coefficient
         if (l == 2){k = stod(s);};
			// the last entry is integer representing the bool
			// for the infinite sink assumption. 1 = true, 0 = false
			if (l == 3){
				integer = stoi(s);
				if (integer == 0){
					infsink = false;
				}
				else {
					infsink = true;
				}
			}
         l++;
      }
		liqIDs.push_back(liqID); surfIDs.push_back(surfID);
		mCoeffs.push_back(k); infSinks.push_back(infsink);

   }
	setWallDeposition(mCoeffs, liqIDs, surfIDs, infSinks);
}

//*****************************************************************************
// Sets a boundary condition in a cell
//
//	@param Loc		Location 
// @param specID  Species ID
// @param bc		BC value [kg/m^3] or [kg/m^2]
//*****************************************************************************
void speciesDriver::setBoundaryCondition(std::string BCType, std::string loc, 
	int specID, double bc){

	int locID = -1;

	if (loc == "north") {locID = 0;};
	if (loc == "south") {locID = 1;};
	if (loc == "east") {locID = 2;};
	if (loc == "west") {locID = 3;};
	assert(locID != -1);

	if (BCType == "dirichlet") {
		setGeneralBoundaryCondition(BCType, locID, specID, bc);
		dummySpec = 1;
	}
	else if (BCType == "newmann"){
		setGeneralBoundaryCondition(BCType, locID, specID, bc);
		dummySpec = 1;
	} 
	else if (BCType == "periodic"){
		setPeriodicBoundaryCondition(locID);
	} 
	else{
		std::string errorMessage = 
			" You have selected an invalid boundary condition ";
		libowskiException::runtimeError(errorMessage);
	}
	
}

//*****************************************************************************
// Sets a boundary condition in a cell for a vector of species
//
//	@param Loc		Location 
// @param ids		vector of species ids
// @param bcs		vector of BC values [kg/m^3] or [kg/m^2]
//*****************************************************************************
void speciesDriver::setBoundaryCondition(std::string BCType, std::string loc, 
	std::vector<int> ids, std::vector<double> bcs){
	int specID;
	double bc;

	// Check to make sure that both vectors are of equal size. 
	// The default arg for bcs is an empty vector. If it is empty is does not
	// do this check
	if (not bcs.empty()){
		assert(ids.size() == bcs.size());
	}

	int locID = -1;

	if (loc == "north") {locID = 0;};
	if (loc == "south") {locID = 1;};
	if (loc == "east") {locID = 2;};
	if (loc == "west") {locID = 3;};
	assert(locID != -1);

	if (BCType == "dirichlet") {
		for (int index = 0; index < ids.size(); index++){
			bc = 0.0;
			specID = ids[index];
			if (not bcs.empty()){ bc = bcs[index];};
			setGeneralBoundaryCondition(BCType, locID, specID, bc);
		}
	}
	else if (BCType == "newmann"){
		for (int index = 0; index < ids.size(); index++){
			bc = 0.0;
			specID = ids[index];
			if (not bcs.empty()){ bc = bcs[index];};
			setGeneralBoundaryCondition(BCType, locID, specID, bc);
		}
	} 
	else if (BCType == "periodic"){
		setPeriodicBoundaryCondition(locID);
	} 
	else{
		std::string errorMessage = 
			" You have selected an invalid boundary condition ";
		libowskiException::runtimeError(errorMessage);
	}
	
}

//*****************************************************************************
// Sets a dirichlet or Newmann boundary condition in a cell
//
// @param type		BC type
//	@param LocID	Location ID
// @param specID  Species ID
// @param bc		BC value [kg/m^3] or [kg/m^2]
//*****************************************************************************
void speciesDriver::setGeneralBoundaryCondition(std::string type, int locID, 
	int specID, double bc){
	int xCellMax = modelPtr->numOfxCells - 1;
	int xCellMin = 0;
	int yCellMax = modelPtr->numOfyCells - 1;
	int yCellMin = 0;

	switch(locID){

		// North location
		case 0: {
			for (int i = 0; i < modelPtr->numOfxCells; i++){
				meshCell* cell = modelPtr->getCellByLoc(i,yCellMax);
				connection* northCon = cell->getConnection(0);
				northCon->boundary = true;
				northCon->boundaryType = type;
				surface* northSurface = northCon->getSurface();
   			species* spec = northSurface->getSpeciesPtr(specID);
				spec->bc = bc;
			}
			break;
		}
		// South location
		case 1: {
			for (int i = 0; i < modelPtr->numOfxCells; i++){
				meshCell* cell = modelPtr->getCellByLoc(i,yCellMin);
				connection* southCon = cell->getConnection(1);
				southCon->boundary = true;
				southCon->boundaryType = type;
				surface* southSurface = southCon->getSurface();
   			species* spec = southSurface->getSpeciesPtr(specID);
				spec->bc = bc;
			}
			break;
		}
		// East location
		case 2: {
			for (int j = 0; j < modelPtr->numOfyCells; j++){
				meshCell* cell = modelPtr->getCellByLoc(xCellMax,j);
				connection* eastCon = cell->getConnection(2);
				eastCon->boundary = true;
				eastCon->boundaryType = type;
				surface* eastSurface = eastCon->getSurface();
   			species* spec = eastSurface->getSpeciesPtr(specID);
				spec->bc = bc;
			}
			break;
		}
		// West location
		case 3: {
			for (int j = 0; j < modelPtr->numOfyCells; j++){
				meshCell* cell = modelPtr->getCellByLoc(xCellMin,j);
				connection* westCon = cell->getConnection(3);
				westCon->boundary = true;
				westCon->boundaryType = type;
				surface* westSurface = westCon->getSurface();
   			species* spec = westSurface->getSpeciesPtr(specID);
				spec->bc = bc;
			}
			break;
		}
	}
}


//*****************************************************************************
// Sets a periodic boundary condition in a cell
//
//	@param LocID	Location ID
//*****************************************************************************
void speciesDriver::setPeriodicBoundaryCondition(int locID){
	int xCellMax = modelPtr->numOfxCells - 1;
	int xCellMin = 0;
	int yCellMax = modelPtr->numOfyCells - 1;
	int yCellMin = 0;

	switch(locID){

		// North location
		case 0: {
			for (int i = 0; i < modelPtr->numOfxCells; i++){
				meshCell* thisCell = modelPtr->getCellByLoc(i,yCellMax);
				meshCell* southCell = modelPtr->getCellByLoc(i,yCellMin);
				connection* thisCellCon = thisCell->getConnection(0);
				assert(thisCellCon->loc == 0);
				thisCellCon->connectionCellPtr = southCell;
			}
			break;
		}
		// South location
		case 1: {
			for (int i = 0; i < modelPtr->numOfxCells; i++){
				meshCell* thisCell = modelPtr->getCellByLoc(i,yCellMin);
				meshCell* northCell = modelPtr->getCellByLoc(i,yCellMax);
				connection* thisCellCon = thisCell->getConnection(1);
				assert(thisCellCon->loc == 1);
				thisCellCon->connectionCellPtr = northCell;
			}
			break;
		}
		// East location
		case 2: {
			for (int j = 0; j < modelPtr->numOfyCells; j++){
				meshCell* thisCell = modelPtr->getCellByLoc(xCellMax,j);
				meshCell* westCell = modelPtr->getCellByLoc(xCellMin,j);
				connection* thisCellCon = thisCell->getConnection(2);
				assert(thisCellCon->loc == 2);
				thisCellCon->connectionCellPtr = westCell;
			}
			break;
		}
		// West location
		case 3: {
			for (int j = 0; j < modelPtr->numOfyCells; j++){
				meshCell* thisCell = modelPtr->getCellByLoc(xCellMin,j);
				meshCell* eastCell = modelPtr->getCellByLoc(xCellMax,j);
				connection* thisCellCon = thisCell->getConnection(3);
				assert(thisCellCon->loc == 3);
				thisCellCon->connectionCellPtr = eastCell;
			}
			break;
		}
	}
}
//*****************************************************************************
// Solves the transient species transport equation
//*****************************************************************************
void speciesDriver::solve(double solveTime){
	VectorD sol;
	MatrixD dA;
	bool augmented = true;
	double timeStep = solveTime - lastSolveTime;
	VectorD defSourceOld, defSourceNew;
	std::ofstream outputFile;

	A = buildTransMatrix(augmented, 0.0);
	if (mpi.rank == 0){
		//dA = Eigen::MatrixXd(A);
		//std::cout << dA.rows() << " " << dA.cols() << std::endl;
		//std::cout << dA.eigenvalues() << std::endl;
		//std::ofstream outputFile;
		//outputFile.open("eigenvalues.txt");
		//outputFile << dA.eigenvalues() << std::endl;
		//outputFile.precision(16);
		//outputFile.setf(ios::fixed);
		//outputFile.setf(ios::showpoint);
		//std::cout << " "  << std::endl;
		//std::cout << dA  << std::endl;
		//std::cout << dA.determinant() << std::endl;
		//std::cout << dA.norm() << std::endl;
		//matrixInit = true;
	}
	N0 = buildInitialConditionVector(augmented);
	sol = expSolver->apply(A, N0, timeStep);
	unpackSolution(sol);
	lastSolveTime = solveTime;
	step += 1;
}

//*****************************************************************************
// Solves the implicit transient species transport equation
//*****************************************************************************
void speciesDriver::solveImplicit(double solveTime){
	VectorD b;
	VectorD sol;
	VectorD solOld;
	MatrixD dA;
	bool augmented = true;
	SparseLU<SparseMatrixD, COLAMDOrdering<int> > LinearSolver;
	double timeStep = solveTime - lastSolveTime;
	double rtol = 1.e-5, diff = 5.0;
	VectorD defSourceOld, defSourceNew;

	A = buildTransMatrix(augmented, 0.0);
	//matrixInit = true;
	solOld = buildInitialConditionVector(augmented);
	sol = intSolver->integrate(A, solOld, timeStep);
	unpackSolution(sol);
	lastSolveTime = solveTime;
}

//*****************************************************************************
// Solves the steady state species transport equation
//*****************************************************************************
void speciesDriver::solve(){
	VectorD sol;
	VectorD b;
	MatrixD dA;
	bool augmented = false;
	SparseLU<SparseMatrixD, COLAMDOrdering<int> > LinearSolver;

	if (mpi.rank == 0){
		A = buildTransMatrix(augmented, 0.0);
		b = buildbVector();
		dA = Eigen::MatrixXd(A);

		LinearSolver.compute(A);
		sol = LinearSolver.solve(b);
		unpackSolution(sol);
	}
}

//*****************************************************************************
// Builds the transition matrix 
//
// @param Augmented	If set to true then the source terms will be moved to into
//							the A matrix and a dummy species will be used to how 
//							their coefficient. If set to false then the source terms
//							will be moved to the other side of the equal sign. This
//							is done for a steady state solve.
//
//	@param dt			Time step over the solve. Only != to zero for implicit
//							transient solves
//*****************************************************************************
SparseMatrixD speciesDriver::buildTransMatrix(bool Augmented, double dt){
	// i, j index of transition matrix
	int i, j;
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	int totalSpecs = numOfSpecs;
	int totalCells = modelPtr->numOfTotalCells;
	int nonZeros = totalCells*totalSpecs*totalSpecs;
	double diffusionCoeff = 0.0;
	double rCon, psi, a, tran, thisCoeff, conDirection, conDist;
	double aP, ab, coeff, defCor, thisSpecSource;
	tripletList.reserve(nonZeros);	

	// if the matrix is not augmented then we no longer need to add a dummy
	// spec to hold the constant source coefficients.
	if (not Augmented){
		dummySpec = 0;
	}

	// Init A matrix
	SparseMatrixD A(totalCells*totalSpecs + dummySpec, totalCells*totalSpecs + 
		dummySpec);
	// Loop over cells
	for (int cellID = 0; cellID < totalCells; cellID++){
		// Gets cell pointer
		meshCell* thisCellPtr = modelPtr->getCellByLoc(cellID);

		// Loop over species
		for (int specID = 0; specID < totalSpecs; specID++){
			// Gets the species pointer
			species* thisSpecPtr = thisCellPtr->getSpecies(specID);
			diffusionCoeff = thisSpecPtr->D;
			thisSpecSource = thisSpecPtr->s;
			// Gets the i matrix index
			i = getAi(cellID, totalCells, specID, totalSpecs);
			// Sets the ap coeff
			aP = 0.0;
			// Set coefficients to zero
			tran = 0.0;
			// i,i coefficient
			thisCoeff = 0.0;

			// loop over cell connections
			for (int conCount = 0; conCount < thisCellPtr->connections.size(); conCount ++){
				// Sets the deferred correction source
				defCor = 0.0;
				// Sets the matrix coefficient to zero
				a = 0.0;
				// Sets a variable that holdes if the matrix coeff is used for setting
				// the boundary condition
				ab = 0.0;
				// Gets cell connection object pointer
				connection* thisCon = thisCellPtr->getConnection(conCount);
				// Gets pointer to connected cell
				meshCell* otherCellPtr = thisCon->connectionCellPtr;
				// Gets the direction required to multiply by the convection transition
				conDirection = thisCon->direction;
				// Gets the distance from this cell center to connection cell center
				conDist = thisCon->distance;

				// If the connection distance is zero then there is no cell in that 
				// direction and the transition rate is zero. This will happen 
				// when modeling 1D cases.
				if (conDist and thisSpecPtr->transport){
					// Computes the transition coefficient for convection
					tran = thisCon->connectionFacePtr->vl*thisCon->area/thisCellPtr->volume;
					// matrix coefficient
					a = std::max(conDirection*tran, 0.0) + 
						diffusionCoeff*thisCon->area/thisCellPtr->volume/conDist;
					// sets the flow coefficients if the pointer is not null
					if (otherCellPtr){
						j = getAi(otherCellPtr->absIndex, totalCells, specID, 
							totalSpecs);
						tripletList.push_back(T(i, j, a));
					}
					// Sets the deferred correction for second order flux
					if (thisCellPtr->secondOrderFlux and not thisCon->boundary){
						// set source for deffered correction
						defCor = calcDefCor(thisCellPtr, thisCon, specID, tran);
						thisSpecSource += defCor;
					}
				}
				// Added sthe Dirichlet boundary condition to the sourse term
				if (thisCon->boundaryType == "dirichlet"){
					surface* thisSurface = thisCon->getSurface();
   				species* surfaceSpecPtr = thisSurface->getSpeciesPtr(specID);
					if (tran){
						thisSpecSource += tran*conDirection*surfaceSpecPtr->bc;
					}
					else{
						thisSpecSource += 2.*a*surfaceSpecPtr->bc;
						ab = -a;
					}
				}
				// Added sthe Newmann boundary condition to the sourse term
				else if (thisCon->boundaryType == "newmann"){
					surface* thisSurface = thisCon->getSurface();
   				species* surfaceSpecPtr = thisSurface->getSpeciesPtr(specID);
					thisSpecSource -= a*surfaceSpecPtr->bc*conDist;
					ab = a;
				}

				// aP coefficient
				aP += (-a + ab);
			}
			// Sets the coefficients for linear source terms
			for (int phyModel = 0; phyModel < thisSpecPtr->sourceTerms.size(); phyModel++){
				for (int specCounter = 0; specCounter < totalSpecs; specCounter++){
					coeff = thisSpecPtr->getTransitionCoeff(specCounter, phyModel,
						thisCellPtr->getScalarData());
					if (specCounter == specID){
						thisCoeff += coeff;
					}
					else{
						j = getAi(cellID, totalCells, specCounter, totalSpecs);
						tripletList.push_back(T(i, j, coeff));
					}
				}
			}

			// Steady state or matrix exp solve
			if (dt == 0.0){thisCoeff += aP;};

			tripletList.push_back(T(i, i, thisCoeff));

			// Sets the constant source terms
			if (Augmented){
				tripletList.push_back(T(i, A.cols()-1, thisSpecSource));
			}
		}
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	return A;
}

//*****************************************************************************
// Builds the initial condition vector
//*****************************************************************************
VectorD speciesDriver::buildInitialConditionVector(bool augmented){
	int i;
	int totalSpecs = numOfSpecs;
	int totalCells = modelPtr->numOfTotalCells;
	if (not augmented){dummySpec = 0;};
	VectorD N0(totalSpecs*totalCells + dummySpec);

	// Loops over cells
	for (int cellID = 0; cellID < totalCells; cellID++){
		// Gets cell pointer
		meshCell* thisCellPtr = modelPtr->getCellByLoc(cellID);

		// Loop over species
		for (int specID = 0; specID < totalSpecs; specID++){
			// Gets the species pointer
			species* thisSpecPtr = thisCellPtr->getSpecies(specID);
			i = getAi(cellID, totalCells, specID, totalSpecs);
			N0[i] = thisSpecPtr->c;
		}
	}
	if(dummySpec){N0[N0.size()-1] = 1.0;};
	return N0;
}

//*****************************************************************************
// Builds the b vector 
//*****************************************************************************
VectorD speciesDriver::buildbVector(){
	int i;
	int totalSpecs = numOfSpecs;
	int totalCells = modelPtr->numOfTotalCells;
	VectorD b(totalSpecs*totalCells);

	// Loops over cells
	for (int cellID = 0; cellID < totalCells; cellID++){
		// Gets cell pointer
		meshCell* thisCellPtr = modelPtr->getCellByLoc(cellID);

		// Loop over species
		for (int specID = 0; specID < totalSpecs; specID++){
			// Gets the species pointer
			species* thisSpecPtr = thisCellPtr->getSpecies(specID);
			i = getAi(cellID, totalCells, specID, totalSpecs);
			b[i] = -thisSpecPtr->s;
		}
	}
	return b;
}

//*****************************************************************************
// Unpacks the solution from the matrix exp solve
//*****************************************************************************
void speciesDriver::unpackSolution(const VectorD& sol){
	int i;
	int totalSpecs = numOfSpecs;
	int totalCells = modelPtr->numOfTotalCells;

	// Loops over cells
	for (int cellID = 0; cellID < totalCells; cellID++){
		// Gets cell pointer
		meshCell* thisCellPtr = modelPtr->getCellByLoc(cellID);

		// Loop over species
		for (int specID = 0; specID < totalSpecs; specID++){
			// Gets the species pointer
			species* thisSpecPtr = thisCellPtr->getSpecies(specID);
			i = getAi(cellID, totalCells, specID, totalSpecs);
			double val = sol[i];
			thisSpecPtr->c = val;

		}
	}
}
//*****************************************************************************
// Calculates the TVD deferred correction
//
// @param cellPtr	pointer to cell
// @param cellCon pointer to cell connection
// @param specID	Species ID
// @pram tran		Transton value	[1/s]
//*****************************************************************************
double speciesDriver::calcDefCor(meshCell* cellPtr, connection* cellCon, 
	int specID, double tran){
	int alphal = 1;
	double rohP = 0.0, rohN = -1.0, rohE = 0.0, rohS = 0.0, rohW = 0.0;
	double rohNN = 0.0, rohSS = 0.0, rohEE = 0.0, rohWW = 0.0;
	double rohbc = 0.0, dir = 0.0;
	double r = 0.0, defCor = 0.0, psi = 0.0;
	double eps = 1.e-16;
	meshCell* northCell = nullptr;
	meshCell* southCell = nullptr;
	meshCell* eastCell = nullptr;
	meshCell* westCell = nullptr;
	meshCell* northNorthCell = nullptr;
	meshCell* southSouthCell = nullptr;
	meshCell* eastEastCell = nullptr;
	meshCell* westWestCell = nullptr;
	meshCell* otherCell = nullptr;
	if (tran < 0.0) {alphal = 0;};
	if (tran == 0.0) {return 0.0;};

	rohP = cellPtr->getSpecCon(specID);
	otherCell = cellCon->connectionCellPtr;

	switch(cellCon->loc){

		// North location
		case 0: {
			northCell = cellPtr->getConnection(0)->connectionCellPtr;
			southCell = cellPtr->getConnection(1)->connectionCellPtr;
			rohN = northCell->getSpecCon(specID);
			northNorthCell = northCell->getConnection(0)->connectionCellPtr;
			// Needed for positive flow
			if (southCell){
				rohS = southCell->getSpecCon(specID);
			}
			else{
				if (cellPtr->getConnection(1)->boundary){
					rohbc = cellPtr->getConnection(1)->getSurface()
						->getSpeciesPtr(specID)->bc;
					if (cellPtr->getConnection(1)->boundaryType == "dirichlet"){
						rohS = 2.*rohbc - rohP;
					}
					else if (cellPtr->getConnection(1)->boundaryType == "newmann"){
						rohS = rohP - rohbc*cellPtr->getConnection(1)->distance;
					}
				}	
			}
			// Needed for negative flow
			if (northNorthCell){
				rohNN = northNorthCell->getSpecCon(specID);
			}
			else{
				if (northCell->getConnection(0)->boundary){
					rohbc = northCell->getConnection(0)->getSurface()
						->getSpeciesPtr(specID)->bc;
					if (northCell->getConnection(0)->boundaryType == "dirichlet"){
						rohNN = 2.*rohbc - rohP;
					}
					else if (northCell->getConnection(0)->boundaryType == "newmann"){
						rohNN = rohP - rohbc*cellCon->distance;
					}
				}
			}
			r = alphal*(rohP - rohS)/(rohN - rohP) + 
				(1.-alphal)*(rohNN - rohN)/(rohN - rohP);
			psi = fluxLim.getPsi(r);
			defCor = 0.5*tran*((1.-alphal)*psi - alphal*psi)*(rohN-rohP);
			break;
		}

		// South location
		case 1: {
			northCell = cellPtr->getConnection(0)->connectionCellPtr;
			southCell = cellPtr->getConnection(1)->connectionCellPtr;
			rohS = southCell->getSpecCon(specID);
			southSouthCell = southCell->getConnection(1)->connectionCellPtr;
			// Needed for positive flow
			if (southSouthCell){
				rohSS = southSouthCell->getSpecCon(specID);
			}
			else{
				if (southCell->getConnection(1)->boundary){
					rohbc = southCell->getConnection(1)->getSurface()
						->getSpeciesPtr(specID)->bc;
					if(southCell->getConnection(1)->boundaryType == "dirichlet"){
						rohSS = 2.*rohbc - rohP;
					}
					else if(southCell->getConnection(1)->boundaryType == "newmann"){
						rohSS = rohP - rohbc*southCell->getConnection(1)->distance;
					}
				}
			}
			// Needed for negative flow
			if (northCell){
				rohN = northCell->getSpecCon(specID);		
			}
			else{
				if (cellPtr->getConnection(0)->boundary){
					rohbc = cellPtr->getConnection(0)->getSurface()
						->getSpeciesPtr(specID)->bc;
					if (cellPtr->getConnection(0)->boundaryType == "dirichlet"){
						rohN = 2.*rohbc - rohP;
					}
					else if (cellPtr->getConnection(0)->boundaryType == "newmann"){
						rohN = rohP - rohbc*cellCon->distance;
					}
				}
			}
			r = alphal*(rohS - rohSS)/(rohP - rohS) + 
				(1.-alphal)*(rohN - rohS)/(rohP - rohS);
			psi = fluxLim.getPsi(r);
			defCor = 0.5*tran*(alphal*psi - (1.-alphal)*psi)*(rohP-rohS);
			break;
		}

		// East location
		case 2: {
			eastCell = cellPtr->getConnection(2)->connectionCellPtr;
			westCell = cellPtr->getConnection(3)->connectionCellPtr;
			rohE = eastCell->getSpecCon(specID);
			eastEastCell = eastCell->getConnection(2)->connectionCellPtr;
			// Needed for positive flow
			if(westCell){
				rohW = westCell->getSpecCon(specID);
			}
			else{
				if (cellPtr->getConnection(3)->boundary){
					rohbc = cellPtr->getConnection(3)->getSurface()
						->getSpeciesPtr(specID)->bc;
					if (cellPtr->getConnection(3)->boundaryType == "dirichlet"){
						rohW = (2.*rohbc - rohP);
					}
					else if (cellPtr->getConnection(3)->boundaryType == "newmann"){
						rohW = rohP - rohbc*cellCon->distance;
					}
				}
			}
			// Needed for negative flow
			if(eastEastCell){
				rohEE = eastEastCell->getSpecCon(specID);
			}
			else{
				if (eastCell->getConnection(2)->boundary){
					rohbc = eastCell->getConnection(2)->getSurface()
						->getSpeciesPtr(specID)->bc;
					if (eastCell->getConnection(2)->boundaryType == "dirichlet"){
						rohEE = (2.*rohbc - rohP);
					}
					else if (eastCell->getConnection(2)->boundaryType == "newmann"){
						rohEE = rohP - rohbc*cellCon->distance;
					}
				}
			}
			r = alphal*(rohP - rohW)/(rohE - rohP) + 
				(1.- alphal)*(rohEE - rohE)/(rohE - rohP);
			psi = fluxLim.getPsi(r);
			defCor = 0.5*tran*((1.-alphal)*psi - alphal*psi)*(rohE-rohP);
			break;
		}

		// West location
		case 3: {
			eastCell = cellPtr->getConnection(2)->connectionCellPtr;
			westCell = cellPtr->getConnection(3)->connectionCellPtr;
			rohW = westCell->getSpecCon(specID);
			westWestCell = westCell->getConnection(3)->connectionCellPtr;
			// Needed for positive flow
			if(westWestCell){
				rohWW = westWestCell->getSpecCon(specID);
			}
			else{
				if (westCell->getConnection(3)->boundary){
					rohbc = westCell->getConnection(3)->getSurface()
						->getSpeciesPtr(specID)->bc;
					if (eastCell->getConnection(3)->boundaryType == "dirichlet"){
						rohWW = (2.*rohbc - rohP);
					}
					else if (eastCell->getConnection(3)->boundaryType == "newmann"){
						rohWW = rohP - rohbc*cellCon->distance;
					}
				}
			}
			// Needed for negative flow
			if(eastCell){
				rohE = eastCell->getSpecCon(specID);
			}
			else{
				if (cellPtr->getConnection(2)->boundary){
					rohbc = cellPtr->getConnection(2)->getSurface()
						->getSpeciesPtr(specID)->bc;
					if (cellPtr->getConnection(2)->boundaryType == "dirichlet"){
						rohE = (2.*rohbc - rohP);
					}
					else if (cellPtr->getConnection(2)->boundaryType == "newmann"){
						rohE = rohP - rohbc*cellCon->distance;
					}
				}
			}
			r = alphal*(rohW - rohWW)/(rohP - rohW) + 
				(1. - alphal)*(rohE - rohP)/(rohP - rohW);
			psi = fluxLim.getPsi(r);
			defCor = 0.5*tran*(alphal*psi - (1.-alphal)*psi)*(rohP-rohW);
			break;
		}

	}
	return defCor;
}

//*****************************************************************************
// Builds a vector of the deferred correction source 
//*****************************************************************************
VectorD speciesDriver::calcDefSourceVector(){
	int i;
	double conDist, conDirection, tran, defCor;
	int totalSpecs = numOfSpecs;
	int totalCells = modelPtr->numOfTotalCells;
	VectorD sourceVector(totalSpecs*totalCells*4);

	// Loop over cells
	i = 0;
	for (int cellID = 0; cellID < totalCells; cellID++){
		// Gets cell pointer
		meshCell* thisCellPtr = modelPtr->getCellByLoc(cellID);

		// Loop over species
		for (int specID = 0; specID < totalSpecs; specID++){
			// Gets the species pointer
			species* thisSpecPtr = thisCellPtr->getSpecies(specID);
			// Set coefficients to zero
			tran = 0.0;

			// loop over cell connections
			for (int conCount = 0; conCount < thisCellPtr->connections.size(); conCount ++){
				// Gets cell connection object pointer
				connection* thisCon = thisCellPtr->getConnection(conCount);
				// Gets pointer to connected cell
				meshCell* otherCellPtr = thisCon->connectionCellPtr;
				// Gets the direction required to multiply by the convection transition
				conDirection = thisCon->direction;
				// Gets the distance from this cell center to connection cell center
				conDist = thisCon->distance;

				// If the connection distance is zero then there is no cell in that 
				// direction and the transition rate is zero. This will happen 
				// when modeling 1D cases.
				if (conDist){
					// Computes the transition coefficient for convection
					tran = thisCon->connectionFacePtr->vl*thisCon->area/thisCellPtr->volume;
					// Sets the deferred correction for second order flux
					if (thisCellPtr->secondOrderFlux and not thisCon->boundary){
						// set source for deffered correction
						defCor = calcDefCor(thisCellPtr, thisCon, specID, tran);
						sourceVector[i] = defCor;
						i += 1;
					}
				}
			}
		}
	}
	return sourceVector;
}

//*****************************************************************************
// Tells solver to rebuild the transition matrix before the next solve
//*****************************************************************************
void speciesDriver::resetMatrix(){
	matrixInit = false;
}

//*****************************************************************************
// Cleans species in the model
//*****************************************************************************
void speciesDriver::clean(){
	numOfSpecs = 0;
	matrixInit = false;
	lastSolveTime = 0.0;
	intSolver->clean();
   for (int i = 0; i < modelPtr->numOfxCells; i++){
      for (int j = 0; j < modelPtr->numOfyCells; j++){
         meshCell* cell = modelPtr->getCellByLoc(i,j);
         cell->cleanSpecies();
      }
   }
}

//*****************************************************************************
// Takes the cell absID and specID and returns a unique location (row 
// or column for that combination in a square matrix
//*****************************************************************************
int speciesDriver::getAi(int cellAbsIndex, int totCells, int specID, int 
		totSpecs){
	int rowsBefore = (cellAbsIndex)*totSpecs;
	return rowsBefore + specID;
}
