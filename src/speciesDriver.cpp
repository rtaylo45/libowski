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
//*****************************************************************************
void speciesDriver::setMatrixExpSolver(std::string solverName, bool krylovFlag,
	int krylovDim){
	expSolver = matrixExponentialFactory::getExpSolver(solverName, krylovFlag,
		krylovDim);
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
// @param molarMass  Molar mass of species [lbm/mol]
// @param [initCon]  Initial concentration [lbm/ft^3]
// @param [diffCoef]	Diffusion coefficient [ft^2/s]
// @param name			Name of the species
//*****************************************************************************
int speciesDriver::addSpecies(double molarMass, double initCon,
	double diffCoeff, std::string name){
   for (int i = 0; i < modelPtr->numOfxCells; i++){
      for (int j = 0; j < modelPtr->numOfyCells; j++){
         meshCell* cell = modelPtr->getCellByLoc(i,j);

         cell->addSpecies(molarMass, initCon, diffCoeff, name);
      }
   }
   int specID = numOfSpecs;
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
		name = result.at(0); mm = stod(result.at(1)); initCon = stod(result.at(2)); 
		D = stod(result.at(3));
		// Adds the species
		specIDs.push_back(addSpecies(mm, initCon, D, name));
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
// Sets the species concentration
//
// @param i       x index
// @param j       y index
// @param specID  Species ID
// @param specCon	Concentration [lbm/ft^3]
//*****************************************************************************
void speciesDriver::setSpeciesCon(int i, int j, int specID, double specCon){
   species* spec = getSpeciesPtr(i, j, specID);
   spec->c = specCon;
}

//*****************************************************************************
// Sets the source terms for a species in a cell
//
// @param i					x index
// @param j					y index
// @param specID			Species ID
// @param coeffs			A vector of species source coefficients
//								[1/s]
// @param transcoeffs	A vector of species souce coefficients for transmutation
//								[ft^2]
// @param s					Constant source in cell [lbm/ft^3/s]
//*****************************************************************************
void speciesDriver::setSpeciesSource(int i, int j, int specID, std::vector<double>
      coeffs, double s, std::vector<double> transCoeffs){
   assert(coeffs.size() == numOfSpecs);
	if (not transCoeffs.empty()){ assert(transCoeffs.size() == numOfSpecs);};
   species* spec = getSpeciesPtr(i, j, specID);
   spec->coeffs = coeffs;
	spec->transCoeffs = transCoeffs;
   spec->s = s;
}

//*****************************************************************************
// Internal function that sets the decay source terms for a species in a cell
//
// @param i					x index
// @param j					y index
// @param specID			Species ID
// @param coeffs			A vector of species source coefficients
//								[1/s]
//*****************************************************************************
void speciesDriver::setDecaySource(int i, int j, int specID, std::vector<double>
      coeffs){
   assert(coeffs.size() == numOfSpecs);
   species* spec = getSpeciesPtr(i, j, specID);
   spec->coeffs = coeffs;
}

//*****************************************************************************
// Internal function that sets the trans source terms for a species in a cell
//
// @param i					x index
// @param j					y index
// @param specID			Species ID
// @param coeffs			A vector of species source coefficients
//								[ft^2]
//*****************************************************************************
void speciesDriver::setTransSource(int i, int j, int specID, std::vector<double>
      coeffs){
   assert(coeffs.size() == numOfSpecs);
   species* spec = getSpeciesPtr(i, j, specID);
   spec->transCoeffs = coeffs;
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
   	   std::vector<double> coeffVect;
   	   int l = 0;
			// Loops though each line to split the string
   	   for (std::string s; iss >> s;){
				// specID is the first entry
   	      if (l == 0){specID = stoi(s);};
				// name is the second entry
   	      if (l == 1){name = s;};
				// the rest of the entries are the source terms
   	      if (l > 1){
   	         coeffVect.push_back(stod(s));
   	      }
   	      l++;
   	   }

			// Loop over cells
			for (int i = 0; i < modelPtr->numOfxCells; i++){
				for (int j = 0; j < modelPtr->numOfyCells; j++){
					if (findex == 0){ 
						setDecaySource(i, j, specID, coeffVect);
					}
					else{
						setTransSource(i, j, specID, coeffVect);
					}
				}
			}
   	}
	}
}

//*****************************************************************************
// Sets a boundary condition in a cell
//
//	@param Loc		Location 
// @param specID  Species ID
// @param bc		BC value [lbm/ft^3]
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
	else if (BCType == "free flow"){
		setGeneralBoundaryCondition(BCType, locID, specID, bc);
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
// @param bc		BC value [lbm/ft^3]
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
	double rtol = 1.e-5, diff;
	VectorD defSourceOld, defSourceNew;

	if (mpi.rank == 0){
		//dA = Eigen::MatrixXd(A*timeStep);
		//std::cout << dA.rows() << " " << dA.cols() << std::endl;
		//std::ofstream outputFile;
		//outputFile.open("matrix.out", std::ios_base::app);
		//outputFile << dA << std::endl;
		//std::cout << dA.eigenvalues() << std::endl;
		//std::cout << " "  << std::endl;
		//std::cout << dA  << std::endl;
		//std::cout << dA.determinant() << std::endl;
		//std::cout << dA.norm() << std::endl;
		//std::cout << N0  << std::endl;
		//matrixInit = true;
	}
	A = buildTransMatrix(augmented, 0.0);
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
				if (conDist){
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
			if (thisSpecPtr->coeffs.size()){
				for (int specCounter = 0; specCounter < totalSpecs; specCounter++){
					coeff = thisSpecPtr->coeffs[specCounter];
					if (not thisSpecPtr->transCoeffs.empty()){ 
						coeff += thisSpecPtr->transCoeffs[specCounter]*thisCellPtr->phi;
					}
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
