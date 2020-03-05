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
}

//*****************************************************************************
// Sets the matrix exponential solver
//
// @param solverName	The name of the materix exponential solver type
//*****************************************************************************
void speciesDriver::setMatrixExpSolver(std::string solverName){
	expSolver = matrixExponentialFactory::getExpSolver(solverName);
}

//*****************************************************************************
// Adds a species to the model
//
// @param molarMass  Molar mass of species [lbm/mol]
// @param [initCon]  Initial concentration [lbm/ft^3]
// @param [diffCoef]	Diffusion coefficient [ft^2/s]
//*****************************************************************************
int speciesDriver::addSpecies(double molarMass, double initCon = 0.0, 
	double diffCoeff = 0.0){
   for (int i = 0; i < modelPtr->numOfxCells; i++){
      for (int j = 0; j < modelPtr->numOfyCells; j++){
         meshCell* cell = modelPtr->getCellByLoc(i,j);

         cell->addSpecies(molarMass, initCon, diffCoeff);
      }
   }
   int specID = numOfSpecs;
   numOfSpecs++;
   return specID;
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
// @param i       x index
// @param j       y index
// @param specID  Species ID
// @param coeffs  A vector of species coefficients size of number of species
//                [lbm/s]
// @param s       Constant source in cell [lbm/ft^3/s]
//*****************************************************************************
void speciesDriver::setSpeciesSource(int i, int j, int specID, std::vector<double>
      coeffs, double s = 0.0){
   assert(coeffs.size() == numOfSpecs);
   species* spec = getSpeciesPtr(i, j, specID);
   spec->coeffs = coeffs;
	if (s != 0.0){dummySpec = 1;};
   spec->s = s;
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

	int xCellMax = modelPtr->numOfxCells - 1;
	int xCellMin = 0;
	int yCellMax = modelPtr->numOfyCells - 1;
	int yCellMin = 0;
	int locID = -1;
	dummySpec = 1;

	if (loc == "north") {locID = 0;};
	if (loc == "south") {locID = 1;};
	if (loc == "east") {locID = 2;};
	if (loc == "west") {locID = 3;};
	assert(locID != -1);


	if (BCType == "dirichlet") {
		setDirichletBoundaryCondition(locID, specID, bc);
	}
	else if (BCType == "periodic"){
		setPeriodicBoundaryCondition(locID);
	} 
	
}

//*****************************************************************************
// Sets a dirichlet boundary condition in a cell
//
//	@param LocID	Location ID
// @param specID  Species ID
// @param bc		BC value [lbm/ft^3]
//*****************************************************************************
void speciesDriver::setDirichletBoundaryCondition(int locID, int specID, 
	double bc){
	int xCellMax = modelPtr->numOfxCells - 1;
	int xCellMin = 0;
	int yCellMax = modelPtr->numOfyCells - 1;
	int yCellMin = 0;

	switch(locID){

		// North location
		case 0: {
			for (int i = 0; i < modelPtr->numOfxCells; i++){
				meshCell* cell = modelPtr->getCellByLoc(i,yCellMax);
				cell->boundary = true;
				cell->boundaryLoc = 0;
   			species* spec = getSpeciesPtr(i, yCellMax, specID);
				spec->bc = bc;
			}
			break;
		}
		// South location
		case 1: {
			for (int i = 0; i < modelPtr->numOfxCells; i++){
				meshCell* cell = modelPtr->getCellByLoc(i,yCellMin);
				cell->boundary = true;
				cell->boundaryLoc = 1;
   			species* spec = getSpeciesPtr(i, yCellMin, specID);
				spec->bc = bc;
			}
			break;
		}
		// East location
		case 2: {
			for (int j = 0; j < modelPtr->numOfyCells; j++){
				meshCell* cell = modelPtr->getCellByLoc(xCellMax,j);
				cell->boundary = true;
				cell->boundaryLoc = 2;
   			species* spec = getSpeciesPtr(xCellMax, j, specID);
				spec->bc = bc;
			}
			break;
		}
		// West location
		case 3: {
			for (int j = 0; j < modelPtr->numOfyCells; j++){
				meshCell* cell = modelPtr->getCellByLoc(xCellMin,j);
				cell->boundary = true;
				cell->boundaryLoc = 3;
   			species* spec = getSpeciesPtr(xCellMin, j, specID);
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
				thisCell->northCellPtr = southCell;
			}
			break;
		}
		// South location
		case 1: {
			for (int i = 0; i < modelPtr->numOfxCells; i++){
				meshCell* thisCell = modelPtr->getCellByLoc(i,yCellMin);
				meshCell* northCell = modelPtr->getCellByLoc(i,yCellMax);
				thisCell->southCellPtr = northCell;
			}
			break;
		}
		// East location
		case 2: {
			for (int j = 0; j < modelPtr->numOfyCells; j++){
				meshCell* thisCell = modelPtr->getCellByLoc(xCellMax,j);
				meshCell* westCell = modelPtr->getCellByLoc(xCellMin,j);
				thisCell->eastCellPtr = westCell;
			}
			break;
		}
		// West location
		case 3: {
			for (int j = 0; j < modelPtr->numOfyCells; j++){
				meshCell* thisCell = modelPtr->getCellByLoc(xCellMin,j);
				meshCell* eastCell = modelPtr->getCellByLoc(xCellMax,j);
				thisCell->westCellPtr = eastCell;
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

	if (not matrixInit){
		A = buildTransMatrix(augmented, 0.0);
		//dA = Eigen::MatrixXd(A*solveTime);
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
		matrixInit = true;
		//N0 = buildInitialConditionVector(augmented);
	}
	N0 = buildInitialConditionVector(augmented);

	sol = expSolver->apply(A, N0, timeStep);
	if (mpi.rank == 0){unpackSolution(sol);};
	lastSolveTime = solveTime;
}

//*****************************************************************************
// Solves the implicit transient species transport equation
//*****************************************************************************
void speciesDriver::solveImplicit(double solveTime){
	VectorD b;
	VectorD sol;
	VectorD cOld;
	MatrixD dA;
	bool augmented = false;
	SparseLU<SparseMatrixD, COLAMDOrdering<int> > LinearSolver;
	double timeStep = solveTime - lastSolveTime;

	cOld = buildInitialConditionVector(augmented);
	A = buildTransMatrix(augmented, timeStep);
	b = -cOld/timeStep + buildbVector();

	LinearSolver.compute(A);
	sol = LinearSolver.solve(b);
	if (mpi.rank == 0){unpackSolution(sol);};
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

	A = buildTransMatrix(augmented, 0.0);
	b = buildbVector();

	LinearSolver.compute(A);
	sol = LinearSolver.solve(-b);
	if (mpi.rank == 0){unpackSolution(sol);};
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
	double aSb = 0.0, aNb = 0.0, aWb = 0.0, aEb = 0.0;
	double rN, rS, rE, rW;
	double psiN, psiS, psiE, psiW;
	double aN, aS, aE, aW;
	double coeff;
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

		// Gets pointer to connecting cells
		meshCell* thisCellNorthCellPtr = thisCellPtr->northCellPtr;
		meshCell* thisCellSouthCellPtr = thisCellPtr->southCellPtr;
		meshCell* thisCellEastCellPtr = thisCellPtr->eastCellPtr;
		meshCell* thisCellWestCellPtr = thisCellPtr->westCellPtr;

		// Gets pointer to the cell faces
		meshCellFace* thisCellNorthFacePtr = thisCellPtr->northFacePtr;
		meshCellFace* thisCellSouthFacePtr = thisCellPtr->southFacePtr;
		meshCellFace* thisCellEastFacePtr = thisCellPtr->eastFacePtr;
		meshCellFace* thisCellWestFacePtr = thisCellPtr->westFacePtr;

		// Gets the transition coefficient for species convective flux
		double nTran = thisCellNorthFacePtr->yVl/thisCellPtr->dy;
		double sTran = thisCellSouthFacePtr->yVl/thisCellPtr->dy;
		double eTran = thisCellEastFacePtr->xVl/thisCellPtr->dx;
		double wTran = thisCellWestFacePtr->xVl/thisCellPtr->dx;

		// Calculates the diffusion transition
		double dN = 1./(thisCellPtr->dy*thisCellPtr->dy); 
		double dS = 1./(thisCellPtr->dy*thisCellPtr->dy);
		double dW = 1./(thisCellPtr->dx*thisCellPtr->dx);
		double dE = 1./(thisCellPtr->dx*thisCellPtr->dx);

		// Loop over species
		for (int specID = 0; specID < totalSpecs; specID++){
			// Gets the species pointer
			species* thisSpecPtr = thisCellPtr->getSpecies(specID);
			diffusionCoeff = thisSpecPtr->D;
			// Gets the i matrix index
			i = getAi(cellID, totalCells, specID, totalSpecs);


			// Conecntration slopes
			rN = calcSpecConvectiveSlope(cellID, specID, 0, nTran);
			rS = calcSpecConvectiveSlope(cellID, specID, 1, sTran);
			rE = calcSpecConvectiveSlope(cellID, specID, 2, eTran);
			rW = calcSpecConvectiveSlope(cellID, specID, 3, wTran);

			// flux limiter
			psiN = fluxLim.getPsi(rN);
			psiS = fluxLim.getPsi(rS);
			psiE = fluxLim.getPsi(rE);
			psiW = fluxLim.getPsi(rW);

			// Matrix Coefficient
			aN = std::max(-nTran,0.0) + psiN/2.*(-std::max(nTran,0.0) -
				std::max(-nTran,0.0)) + diffusionCoeff*dN;
			aS = std::max(sTran,0.0) + psiS/2.*(std::max(-sTran,0.0) -
				std::max(sTran, 0.0)) + diffusionCoeff*dS;
			aE = std::max(-eTran,0.0) + psiE/2.*(-std::max(-eTran,0.0) -
				std::max(eTran, 0.0)) + diffusionCoeff*dE;
			aW = std::max(wTran, 0.0) + psiW/2.*(std::max(-wTran, 0.0) -
				std::max(wTran,0.0)) + diffusionCoeff*dW;

			// Sets the north flow coefficient 
			if (thisCellNorthCellPtr){
				j = getAi(thisCellNorthCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				tripletList.push_back(T(i, j, aN));
			}	
			// Sets the south flow coefficient
			if (thisCellSouthCellPtr){
				j = getAi(thisCellSouthCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				tripletList.push_back(T(i, j, std::max(aS,0.0)));
			}
			// Sets the east flow coefficient
			if(thisCellEastCellPtr){
				j = getAi(thisCellEastCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				tripletList.push_back(T(i, j, std::max(aE,0.0)));
			}
			// Sets the west flow coefficient
			if(thisCellWestCellPtr){
				j = getAi(thisCellWestCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				tripletList.push_back(T(i, j, std::max(aW,0.0)));
			}
			// Sets the coefficients for non-constant source terms
			double thisCoeff = 0.0;
			if (thisSpecPtr->coeffs.size()){
				for (int specCounter = 0; specCounter < totalSpecs; specCounter++){
					coeff = thisSpecPtr->coeffs[specCounter];
					if (specCounter == specID){
						thisCoeff += coeff;
					}
					else{
						j = getAi(cellID, totalCells, specCounter, totalSpecs);
						tripletList.push_back(T(i, j, coeff));
					}
				}
			}

			if (thisCellPtr->boundaryLoc == 0){thisSpecPtr->s += 2.*aN*thisSpecPtr->bc;};
			if (thisCellPtr->boundaryLoc == 1){thisSpecPtr->s += 2.*aS*thisSpecPtr->bc;};
			if (thisCellPtr->boundaryLoc == 2){thisSpecPtr->s += 2.*aE*thisSpecPtr->bc;};
			if (thisCellPtr->boundaryLoc == 3){thisSpecPtr->s += 2.*aW*thisSpecPtr->bc;};

			if (thisCellPtr->boundaryLoc == 0){aNb = aN;};
			if (thisCellPtr->boundaryLoc == 1){aSb = aS;};
			if (thisCellPtr->boundaryLoc == 2){aEb = aE;};
			if (thisCellPtr->boundaryLoc == 3){aWb = aW;};

			// Calculates the x portion for the cell coefficient
			double aPx = -std::max(eTran,0.0) - std::max(-eTran,0.0) +
				psiE/2.*(std::max(eTran,0.0) + std::max(-eTran,0.0)) +
				psiW/2.*(std::max(wTran,0.0) + std::max(-wTran,0.0)) -
				diffusionCoeff*dE - diffusionCoeff*dW;

			// Calculates the y portion for the cell coefficient
			double aPy = -std::max(nTran,0.0) - std::max(-sTran,0.0) + 
				psiN/2.*(std::max(nTran,0.0) + std::max(-nTran,0.0)) +
				psiS/2.*(std::max(sTran,0.0) + std::max(-sTran,0.0)) - 
				diffusionCoeff*dS - diffusionCoeff*dN;

			// Adds the coeff for this species implicit solve
			if (dt != 0.0){thisCoeff += aPx + aPy - 1./dt;};
			// Steady state or matrix exp solve
			if (dt == 0.0){thisCoeff += aPx + aPy;};
			// Adds the coefficents if the cell in a boundary
			thisCoeff -= (aSb + aNb + aWb + aEb);
			tripletList.push_back(T(i, i, thisCoeff));

			// Sets the constant source terms
			if (Augmented){tripletList.push_back(T(i, A.cols()-1, thisSpecPtr->s));};
			// Resets the boundary coefficients
			aSb = 0.0, aNb = 0.0, aWb = 0.0, aEb = 0.0;
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
// Claculatest the convective species slop in a cell
//
// @param cellID	Global cellID
// @param specID	Species ID
// @pram tran		Transton value	[1/s]
// @param loc		Location of adjacent cell
//			0 = north
//			1 = south
//			2 = east
//			3 = west
//*****************************************************************************
double speciesDriver::calcSpecConvectiveSlope(int cellID, int specID, 
		int loc, double tran){
		double alphal = 0.0;
		double rohP = 0.0, rohN = 0.0, rohE = 0.0, rohS = 0.0, rohW = 0.0;
		double rohNN, rohSS, rohEE, rohWW;
		double r;
		if (tran < 0.0) {alphal = 1.0;};
		if (tran == 0.0) {return 0.0;};

		// This cell pointer
		meshCell* thisCellPtr = modelPtr->getCellByLoc(cellID);
	
		// Gets pointer to connecting cells
		meshCell* thisCellNorthCellPtr = thisCellPtr->northCellPtr;
		meshCell* thisCellSouthCellPtr = thisCellPtr->southCellPtr;
		meshCell* thisCellEastCellPtr = thisCellPtr->eastCellPtr;
		meshCell* thisCellWestCellPtr = thisCellPtr->westCellPtr;

		// Gets pointer to conncetions of connections
		meshCell* thisCellNorthNorthCellPtr = nullptr;
		meshCell* thisCellSouthSouthCellPtr = nullptr;
		meshCell* thisCellEastEastCellPtr = nullptr;
		meshCell* thisCellWestWestCellPtr = nullptr;

		if (thisCellNorthCellPtr){thisCellNorthNorthCellPtr = thisCellNorthCellPtr->northCellPtr;};
		if (thisCellSouthCellPtr){thisCellSouthSouthCellPtr = thisCellSouthCellPtr->southCellPtr;};
		if (thisCellEastCellPtr){thisCellEastEastCellPtr = thisCellEastCellPtr->eastCellPtr;};
		if (thisCellWestCellPtr){thisCellWestWestCellPtr = thisCellWestCellPtr->westCellPtr;};

		if (thisCellNorthCellPtr){rohN = thisCellNorthCellPtr->getSpecCon(specID);};
		if (thisCellSouthCellPtr){rohS = thisCellSouthCellPtr->getSpecCon(specID);};
		if (thisCellEastCellPtr){rohE = thisCellEastCellPtr->getSpecCon(specID);};
		if (thisCellWestCellPtr){rohW = thisCellWestCellPtr->getSpecCon(specID);};

		// Gets this cells species concentration
		rohP = thisCellPtr->getSpecCon(specID);

		switch(loc){

			// North location
			case 0: {
				rohNN = (thisCellNorthNorthCellPtr) ? 
					thisCellNorthNorthCellPtr->getSpecCon(specID) : 0.0;
				r = (1.-alphal)*(rohP - rohS)/(rohN - rohP) + 
					alphal*(rohNN - rohN)/(rohN - rohP);
				
				break;
			}

			// South location
			case 1: {
				rohSS = (thisCellSouthSouthCellPtr) ? 
					thisCellSouthSouthCellPtr->getSpecCon(specID) : 0.0;
				r = (1.-alphal)*(rohS - rohSS)/(rohP - rohS) + 
					alphal*(rohN - rohS)/(rohP - rohS);
				break;
			}

			// East location
			case 2: {
				rohEE = (thisCellEastEastCellPtr) ? 
					thisCellEastEastCellPtr->getSpecCon(specID) : 0.0;
				r = (1.-alphal)*(rohP - rohW)/(rohE - rohP) + 
					alphal*(rohEE - rohE)/(rohE - rohP);
				break;
			}

			// West location
			case 3: {
				rohWW = (thisCellWestWestCellPtr) ? 
					thisCellWestWestCellPtr->getSpecCon(specID) : 0.0;
				r = (1.-alphal)*(rohW - rohWW)/(rohP - rohW) + 
					alphal*(rohE - rohP)/(rohP - rohW);
				break;
			}

		}
		return r;
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
