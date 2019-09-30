//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include "speciesDriver.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

//*****************************************************************************
// Constructor
//
// @param model A pointer to the mesh model
//*****************************************************************************
speciesDriver::speciesDriver(modelMesh* model){
	modelPtr = model;
}

//*****************************************************************************
// Adds a species to the model
//
// @param molarMass  Molar mass of species [lbm/mol]
// @param [initCon]  Initial concentration [lbm/ft^3]
//*****************************************************************************
int speciesDriver::addSpecies(double molarMass, double initCon = 0.0){
   for (int i = 0; i < modelPtr->numOfxCells; i++){
      for (int j = 0; j < modelPtr->numOfyCells; j++){
         meshCell* cell = modelPtr->getCellByLoc(i,j);

         cell->addSpecies(molarMass, initCon);
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
// @param i       x index
// @param j       y index
// @param specID  Species ID
// @param bc			BC value [lbm/ft^3]
//*****************************************************************************
void speciesDriver::setBoundaryCondition(int i, int j, int specID, double bc){
   meshCell* cell = modelPtr->getCellByLoc(i,j);
	cell->boundary = true;
	dummySpec = 1;
   species* spec = getSpeciesPtr(i, j, specID);
	spec->bc = bc;
	
}

//*****************************************************************************
// Solves the species transport equation
//*****************************************************************************
void speciesDriver::solve(double solveTime){
	Eigen::VectorXd sol;
	SolverType ExpSolver;
	Eigen::MatrixXd dA;
	double timeStep = solveTime - lastSolveTime;

	if (not matrixInit){
		A = buildTransMatrix();
		//dA = Eigen::MatrixXd(A);
		//std::cout << A  << std::endl;
		//std::cout << dA.eigenvalues() << std::endl;
		//std::cout << dA.determinant() << std::endl;
		//std::cout << N0  << std::endl;
		matrixInit = true;
	}
	N0 = buildInitialConditionVector();

	sol = ExpSolver.solve(A, N0, timeStep);
	unpackSolution(sol);
	lastSolveTime = solveTime;
}

//*****************************************************************************
// Builds the transition matrix 
//*****************************************************************************
Eigen::SparseMatrix<double> speciesDriver::buildTransMatrix(){
	// i, j index of transition matrix
	int i, j;
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	int totalSpecs = numOfSpecs;
	int totalCells = modelPtr->numOfTotalCells;
	int nonZeros = totalCells*totalSpecs*totalSpecs;
	tripletList.reserve(nonZeros);	

	// Init A matrix
	Eigen::SparseMatrix<double> A(totalCells*totalSpecs + dummySpec, 
			totalCells*totalSpecs + dummySpec);
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
		double dN = 0.0; 
		double dS = 0.0;
		double dW = 0.0;
		double dE = 0.0;

		// Loop over species
		for (int specID = 0; specID < totalSpecs; specID++){
			// Gets the species pointer
			species* thisSpecPtr = thisCellPtr->getSpecies(specID);
			// Gets the i matrix index
			i = getAi(cellID, totalCells, specID, totalSpecs);

			// Sets the north flow coefficient 
			if (thisCellNorthCellPtr){
				j = getAi(thisCellNorthCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				double r = calcSpecConvectiveSlope(cellID, specID, 0, nTran);
				double psi = fluxLim.getPsi(r);
				double aN = std::max(-nTran,0.0) + psi/2.*(-std::max(nTran,0.0) -
					std::max(-nTran,0.0)) + dN;
				tripletList.push_back(T(i, j, aN));
			}	
			// Sets the south flow coefficient
			if (thisCellSouthCellPtr){
				j = getAi(thisCellSouthCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				double r = calcSpecConvectiveSlope(cellID, specID, 1, sTran);
				double psi = fluxLim.getPsi(r);
				double aS = std::max(sTran,0.0) + psi/2.*(std::max(-sTran,0.0) -
					std::max(sTran, 0.0)) + dS;
				tripletList.push_back(T(i, j, std::max(aS,0.0)));
			}
			// Sets the east flow coefficient
			if(thisCellEastCellPtr){
				j = getAi(thisCellEastCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				double r = calcSpecConvectiveSlope(cellID, specID, 2, eTran);
				double psi = fluxLim.getPsi(r);
				double aE = std::max(-eTran,0.0) + psi/2.*(-std::max(-eTran,0.0) -
					std::max(eTran, 0.0)) + dE;
				tripletList.push_back(T(i, j, std::max(aE,0.0)));
			}
			// Sets the west flow coefficient
			if(thisCellWestCellPtr){
				j = getAi(thisCellWestCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				double r = calcSpecConvectiveSlope(cellID, specID, 3, wTran);
				double psi = fluxLim.getPsi(r);
				double aW = std::max(wTran, 0.0) + psi/2.*(std::max(-wTran, 0.0) -
					std::max(wTran,0.0)) + dW;
				tripletList.push_back(T(i, j, std::max(aW,0.0)));
			}
			// Sets the coefficients for non-constant source terms
			double thisCoeff = 0.0;
			if (thisSpecPtr->coeffs.size()){
				for (int specCounter = 0; specCounter < totalSpecs; specCounter++){
					double coeff = thisSpecPtr->coeffs[specCounter];
					if (specCounter == specID){
						thisCoeff += coeff;
					}
					else{
						j = getAi(cellID, totalCells, specCounter, totalSpecs);
						tripletList.push_back(T(i, j, coeff));
					}
				}
			}

			double rS = calcSpecConvectiveSlope(cellID, specID, 1, sTran);
			double aSb = std::max(sTran,0.0) + fluxLim.getPsi(rS)*(std::max(-sTran,0.0) -
				std::max(sTran, 0.0)) + dS;
			if (thisCellPtr->boundary){thisSpecPtr->s += 2.*aSb*thisSpecPtr->bc;};
			// Sets the constant source terms
			tripletList.push_back(T(i, A.cols()-1, thisSpecPtr->s));

			// Calculates the x portion for the cell coefficient
			double rW = calcSpecConvectiveSlope(cellID, specID, 3, wTran);
			double rE = calcSpecConvectiveSlope(cellID, specID, 2, eTran);
			double psiW = fluxLim.getPsi(rW);
			double psiE = fluxLim.getPsi(rE);
			double aPx = -std::max(eTran,0.0) - std::max(-eTran,0.0) +
				psiE/2.*(std::max(eTran,0.0) + std::max(-eTran,0.0)) +
				psiW/2.*(std::max(wTran,0.0) + std::max(-wTran,0.0)) -
				dE - dW;

			// Calculates the y portion for the cell coefficient
			double rN = calcSpecConvectiveSlope(cellID, specID, 0, nTran);
			//double rS = calcSpecConvectiveSlope(cellID, specID, 1, sTran);
			double psiN = fluxLim.getPsi(rN);
			double psiS = fluxLim.getPsi(rS);
			double aPy = -std::max(nTran,0.0) - std::max(-sTran,0.0) + 
				psiN/2.*(std::max(nTran,0.0) + std::max(-nTran,0.0)) +
				psiS/2.*(std::max(sTran,0.0) + std::max(-sTran,0.0)) - 
				dS - dN;
			if (thisCellPtr->boundary){aPy -= aSb;};

			// Adds the coeff for this species 
			thisCoeff += aPx + aPy;
			//std::cout << std::max(nTran,0.0)<< " "<< std::max(sTran,0.0) <<std::endl;
			tripletList.push_back(T(i, i, thisCoeff));
		}
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	return A;
}

//*****************************************************************************
// Builds the initial condition vector
//*****************************************************************************
Eigen::VectorXd speciesDriver::buildInitialConditionVector(){
	int i;
	int totalSpecs = numOfSpecs;
	int totalCells = modelPtr->numOfTotalCells;
	Eigen::VectorXd N0(totalSpecs*totalCells + dummySpec);

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
// Unpacks the solution from the matrix exp solve
//*****************************************************************************
void speciesDriver::unpackSolution(Eigen::VectorXd sol){
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
