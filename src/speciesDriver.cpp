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
	cell->solved = true;
   species* spec = getSpeciesPtr(i, j, specID);
	spec->c = bc;
	
}

//*****************************************************************************
// Solves the species transport equation
//*****************************************************************************
void speciesDriver::solve(double solveTime){
	Eigen::VectorXd sol;
	SolverType ExpSolver;

	Eigen::SparseMatrix<double> A = buildTransMatrix();
	Eigen::VectorXd N0 = buildInitialConditionVector();

	sol = ExpSolver.solve(A, N0, solveTime);
	//std::cout << sol;
	unpackSolution(sol);
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
		if (thisCellPtr->solved) continue;

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
		double nTran = -thisCellNorthFacePtr->yVl/thisCellPtr->dy;
		double sTran = thisCellSouthFacePtr->yVl/thisCellPtr->dy;
		double eTran = -thisCellEastFacePtr->xVl/thisCellPtr->dx;
		double wTran = thisCellWestFacePtr->xVl/thisCellPtr->dx;

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
				tripletList.push_back(T(i, j, std::max(nTran,0.0)));
			}	
			// Sets the south flow coefficient
			if (thisCellSouthCellPtr){
				j = getAi(thisCellSouthCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				tripletList.push_back(T(i, j, std::max(sTran,0.0)));
			}
			// Sets the east flow coefficient
			if(thisCellEastCellPtr){
				j = getAi(thisCellEastCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				tripletList.push_back(T(i, j, std::max(eTran,0.0)));
			}
			// Sets the west flow coefficient
			if(thisCellWestCellPtr){
				j = getAi(thisCellWestCellPtr->absIndex, totalCells, specID, 
					totalSpecs);
				tripletList.push_back(T(i, j, std::max(wTran,0.0)));
			}
			// Sets the coefficients for non-constant source terms
			double thisCoeff = 0.0;
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
			// Sets the constant source terms
			tripletList.push_back(T(i, A.cols()-1, thisSpecPtr->s));

			// Adds the coeff for this species 
			thisCoeff += std::min(nTran,0.0) + std::min(sTran,0.0) 
				+ std::min(wTran,0.0) + std::min(eTran,0.0);
			//std::cout << std::max(nTran,0.0)<< " "<< std::max(sTran,0.0) <<std::endl;
			tripletList.push_back(T(i, i, thisCoeff));
		}
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	//std::cout << A << std::endl;
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
	N0[N0.size()-1] = 1.0;
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
