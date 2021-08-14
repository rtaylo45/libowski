#include "mpiProcess.h"
#include "modelMesh.h"
#include "speciesDriver.h"
#include "meshCell.h"
#include "species.h"

//**************************************************************************
// 2D Neutron precuror problem
//
//	Problem equations:
//	dCi/dt = -v*dCi/dx - lambda_i*Ci + beta_i*phi*sin(pi*x/100)*sin(pi*y/50)
//
//	Domaine:
//		x = [0, 4]		m
//		y = [0, 0.5]   m
//		t = [0, 60]		s
//
//		i | lambda_i | beta_i
//		---------------------
//		1 | 0.0127	 | 0.0006 
//   	2 | 0.0317   | 0.00364 
//   	3 | 0.115    | 0.00349 
//   	4 | 0.311    | 0.00628 
//   	5 | 1.4      | 0.00179
//   	6 | 3.87     | 0.0007
//
//	Initial conditions and BC's:
//			C_{1}	= 0.0 
//			C_{2}	= 0.0 
//			C_{3}	= 0.0 
//			C_{4}	= 0.0 
//			C_{5}	= 0.0 
//			C_{6}	= 0.0 
//
//			C_{i}(0, y, t)			= C_{i}(400, y, t)
//			dC_{i}(x, 0, t)/dt	= 0
//			dC_{i}(x, 50, t)/dt	= 0
//
//	Solution:
//		Calculated by matlab
// 
// Note:
//		The xCells and yCells were choosen so that the yCells do not overlap
//		outside of the core region. i.e. each cell is fully in the core or 
//		not in the core
//**************************************************************************
void neutronPrecursors(int myid, std::string solverType){
	double t = 0.0;
	int steps = 1;
	double totalTime = 60.0;
	double dt = totalTime/steps;
	int xCells = 5, yCells = 20;
	double xLength = 50.0, yLength = 400.0;
	double v_y = 100.0, v_x = 0.0;
	double xc, yc, s, g;
	MatrixD refSolData;
	std::vector<int> ids;
	std::vector<double> decay = {0.0127, 0.0317, 0.115, 0.311, 1.4, 3.87};
	std::vector<double> beta = {0.06, 0.364, 0.349, 0.628, 0.179, 0.07};
	std::vector<double> bcs = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	std::string path = getDataPath() + "caseStudy/";
	std::string solutionPath = path + "neutronPrecursors/dt60/";
	std::string speciesNamesFile = path + "neutronPrecursorsInputNames.dat";
	std::string limiter = "First order upwind";
	std::string outputFileName = "caseStudyNeutronPrecursors.out";
	std::string solutionFileName = solutionPath + "solutionVel100.csv";
	std::string outputFileNameMatrix = "caseStudyNeutronprecursors"
		//+solverType+"Substeps4.csv";
		+solverType+".csv";
	FILE * pOutputFile;
	meshCell* cell = nullptr;
	pOutputFile = fopen(outputFileName.c_str(), "a");

	if (myid == 0){
		fprintf(pOutputFile, "Solver: %s \n", solverType.c_str());
		fprintf(pOutputFile, "Total problem time: %4.3e\n", totalTime);
		fprintf(pOutputFile, "yLength: %4.3e\n", yLength);
		fprintf(pOutputFile, "xLength: %4.3e\n", xLength);
		fprintf(pOutputFile, "dt: %8.7e\n", dt);
		fprintf(pOutputFile, "yCells: %d\n", yCells);
		fprintf(pOutputFile, "xCells: %d\n", xCells);
		fprintf(pOutputFile, "%s %s %s %s %s %s %s %s %s \n", "variables", 
			"x", "y", "G1", "G2", "G3", "G4", "G5", "G6");
	}
	fclose(pOutputFile);

	// Builds the model mesh
	modelMesh model(xCells, yCells, xLength, yLength);

	// set velocity
	model.setConstantYVelocity(v_y);
	model.setConstantXVelocity(v_x);

	// Adds boundary surface
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");

	// Inits the species driver
	speciesDriver spec = speciesDriver(&model);

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);
	//spec.setIntegratorSolver("explicit", solverType);

	// Sets the flux limiter type
	spec.setFluxLimiter(limiter);

	// Species IDs
	ids = spec.addSpeciesFromFile(speciesNamesFile);

	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	readCSV(refSolData, solutionFileName);

	// Adds the boundary conditions
	spec.setBoundaryCondition("newmann","east", ids, bcs);
	spec.setBoundaryCondition("newmann","west", ids, bcs);
	spec.setBoundaryCondition("periodic","south", ids);
	spec.setBoundaryCondition("periodic","north", ids);

	// Sets the source terms
	for (int id = 0; id < ids.size(); id++){
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				std::vector<double> coeffs(ids.size(), 0.0);
				coeffs[id] = -decay[id];
				meshCell* cell = model.getCellByLoc(i,j);
				double y = cell->y, x = cell->x;
				double y1 = y - model.dy/2., x1 = x - model.dx/2.;
				double y2 = y + model.dy/2., x2 = x + model.dx/2.;
				double sy = (1./model.dy)*(100./M_PI)*(cos(M_PI*y1/100.) - 
				cos(M_PI*y2/100.));
				double sx = (1./model.dx)*(50./M_PI)*(cos(M_PI*x1/50.) - 
				cos(M_PI*x2/50.));

				if (y < 100.){
					s = 1.e13*sy*sx;
				}
				else{
					s = 0.0;
				}
				model.setCellNeutronFlux(i, j, s);
				spec.setSpeciesSource(i, j, id, coeffs, beta[id]*s);
				//std::cout << i << " " << j << " " << id << " " << s << std::endl;
			}
		}
	}


	if (myid == 0){
		spec.writeTransitionMatrixToFile("transitionMatrixNeutronPrecursors.csv");
		spec.writeInitialConditionToFile("initialConditionNeutronPrecursors.csv");
	}


	pOutputFile = fopen(outputFileName.c_str(), "a");
	// Loops to solve the problem
	double totalSolvetime = 0.0;
	for (int k = 0; k < steps; k++){
		t = t + dt;
		auto start = std::chrono::high_resolution_clock::now();
		spec.solve(t);
		//spec.solveImplicit(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		//std::cout << "Time: " << t << " Solve Time: "
		//<< duration.count()/1.e6 << std::endl;
		if (myid == 0){
			totalSolvetime += duration.count()/1.e6;
		}	
	//}
		if (myid == 0){
			int index = 0;
			// Loops to print results
			fprintf(pOutputFile, "time: %5.4e\n", t);
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j);
					xc = cell->x;
					yc = cell->y;
					fprintf(pOutputFile, "%5.4e %5.4e ", xc, yc);
					// Loops over the species 
					for (int id = 0; id < ids.size(); id++){
						double con = spec.getSpecies(i, j, ids[id]);
						fprintf(pOutputFile, "%17.16e ", con);
						solData(index, k) = con;
						index ++;
					}
					fprintf(pOutputFile, "\n");
				}
			}
			fprintf(pOutputFile, "\n");
			VectorD refSol = refSolData.col(k);
			VectorD sol = solData.col(k);
			double rmse = computeRelativeRMSE(refSol, sol);
			double E1 = computeRelativeE1(refSol, sol);
			double E2 = computeRelativeE2(refSol, sol);
			double Einf = computeRelativeEinfty(refSol, sol);
			printf("%s %f %e %e %e\n", solverType.c_str(), t, Einf, E1, E2);
			//printf("%e %e %e", Einf, E1, E2);
			//if (solverType != "pade-method1"){
			//	assert(rmse < 1.e-10);
			//}
			//else{
			//	assert(rmse < 1.e-2);
			//}
		}
	}
	if (myid == 0){
		fprintf(pOutputFile, "end\n");
		//fclose(pOutputFile);
		std::cout << solverType << " Solve time: " << totalSolvetime << std::endl;
		writeCSV(solData, outputFileNameMatrix);
	}
}
//*****************************************************************************
// MSR lump depletion
//*****************************************************************************
void msrLumpDepletion(int myid, std::string solverType){
	double t = 0.0;
	int steps = 20;
	double depletionTime = 200.; // Days
	//double depletionTime = 20.; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 1, yCells = 1;
	double xLength = 1., yLength = 1., xc, yc;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string solPath = getDataPath() + "caseStudy/msrLumpDepletion/";
	std::string outputFileName = "msrMicroLumpDepletion.out";
	std::string outputFileNameMatrix = solverType+"MSRLumpDepletionMicro.csv";
	std::string limiter = "First order upwind";
	std::ofstream outputFile;
	outputFile.open(outputFileName, std::ios_base::app);
	MatrixD refSolData;
	FILE * pOutputFile;
	meshCell* cell = nullptr;
	pOutputFile = fopen(outputFileName.c_str(), "a");

	if (myid == 0){
		fprintf(pOutputFile, "Solver: %s \n", solverType.c_str());
		fprintf(pOutputFile, "Total problem time: %4.3e\n", totalTime);
		fprintf(pOutputFile, "yLength: %4.3e\n", yLength);
		fprintf(pOutputFile, "xLength: %4.3e\n", xLength);
		fprintf(pOutputFile, "dt: %8.7e\n", dt);
		fprintf(pOutputFile, "yCells: %d\n", yCells);
		fprintf(pOutputFile, "xCells: %d\n", xCells);
		fprintf(pOutputFile, "%s %s %s %s \n", "variables", 
			"x", "y", "Name");
	fclose(pOutputFile);
	}

	// Sets file paths for the input data	
	std::string speciesNamesFile = path + "speciesInputNamesMSRMicro.dat";
	std::string speciesDecayFile = path + "speciesInputDecayMSRMicro.dat";
	std::string speciesTransFile = path + "speciesInputTransMSRMicro.dat";

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

	// Writes transition matrix and initial condition
	spec.writeTransitionMatrixToFile("transitionMatrixMSRMicroLumpDepletion.csv");
	spec.writeInitialConditionToFile("initialConditionMSRMicroLumpDepletion.csv");
	pOutputFile = fopen(outputFileName.c_str(), "a");

	// Gets the solution
	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	//readCSV(refSolData, std::string(solPath + "solutionMSRSmallLumpDepletion.csv"));

	// Loops over the time steps to solve
	double totalSolveTime = 0;
	for (int k = 0; k < steps; k++){
		auto start = std::chrono::high_resolution_clock::now();
		t = t + dt;
		spec.solve(t);
		//spec.solveImplicit(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		solveTime = duration.count()/1.e6;
		totalSolveTime += solveTime;
		if (myid == 0){

			cell = model.getCellByLoc(0,0); xc = cell->x; yc = cell->y;

			fprintf(pOutputFile, "time: %5.4e\n", t);
			//printf("time: %5.4e\n", t);

			for (int id = 0; id < ids.size(); id++){
				std::string name = spec.getSpeciesName(ids[id]);
				double con = spec.getSpecies(0, 0, ids[id]);
				fprintf(pOutputFile, "%5.4e %5.4e %s %17.16e \n", xc, yc, name.c_str(), con);
				printf("%5.4e %5.4e %s %17.16e \n", xc, yc, name.c_str(), con);
				solData(id, k) = con;
			}
			fprintf(pOutputFile, "\n");
			//VectorD refSol = refSolData.col(k), sol = solData.col(k);
			//double E1 = computeRelativeE1(refSol, sol);
			//double E2 = computeRelativeE2(refSol, sol);
			//double Einf = computeRelativeEinfty(refSol, sol);
			//printf("%s %e %e %e %e\n", solverType.c_str(), t, Einf, E1, E2);
			//printf("Solve Step: %d %s Solve Time: %f RMSE %e\n", k, solverType.c_str(), solveTime, rmse);
		}
	}
	if (myid == 0){
		fprintf(pOutputFile, "end\n");
		//double rmse = computeRelativeRMSE(refSolData, solData);
		//double rmse = 0.0;
		printf("%s Solve Time: %f \n", solverType.c_str(), totalSolveTime);
		//writeCSV(solData, outputFileNameMatrix);
	}
}

//*****************************************************************************
// MSR 2-D depletion 
//
// No mass transport this is for the M&C paper
//*****************************************************************************
void msr2DDepletion(int myid, std::string solverType){
	double t = 0.0;
	int steps = 10;
	double coreLength = 1.92024;
	double xLength = 0.6858, yLength = 3.*coreLength;		// Meters
	double v_y = 0.25;											// m/s
	double depletionTime = 200.;							// Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 9, yCells = 27;
	double xc, yc, s;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string outputFileName = "msr2DDepletionSmall.out";
	std::string limiter = "First order upwind";
	std::string outputFileNameMatrix = solverType+"MSR2DDepletion.csv";
	std::ofstream outputFile;
	outputFile.open(outputFileName, std::ios_base::app);
	MatrixD refSolData;
	FILE * pOutputFile;
	meshCell* cell = nullptr;

	if (myid == 0){
		pOutputFile = fopen(outputFileName.c_str(), "a");
		fprintf(pOutputFile, "Solver: %s \n", solverType.c_str());
		fprintf(pOutputFile, "Total problem time: %4.3e\n", totalTime);
		fprintf(pOutputFile, "yLength: %4.3e\n", yLength);
		fprintf(pOutputFile, "xLength: %4.3e\n", xLength);
		fprintf(pOutputFile, "dt: %8.7e\n", dt);
		fprintf(pOutputFile, "yCells: %d\n", yCells);
		fprintf(pOutputFile, "xCells: %d\n", xCells);
		fclose(pOutputFile);
	}

	// Sets file paths for the input data	
	std::string speciesNamesFile = path + "speciesInputNamesMSRSmall.dat";
	std::string speciesDecayFile = path + "speciesInputDecayMSRSmall.dat";
	std::string speciesTransFile = path + "speciesInputTransMSRSmall.dat";

	// Builds the model
	modelMesh model(xCells, yCells, xLength, yLength);

	// Adds boundary surface
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");

	// set y velocity
	model.setConstantYVelocity(v_y);

	// Builds the species object
	speciesDriver spec = speciesDriver(&model);

	// Sets the source terms
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			std::vector<double> coeffs(ids.size(), 0.0);
			meshCell* cell = model.getCellByLoc(i,j);
			double y = cell->y, x = cell->x;
			double y1 = y - model.dy/2., x1 = x - model.dx/2.;
			double y2 = y + model.dy/2., x2 = x + model.dx/2.;
			double sy = (1./model.dy)*(coreLength/M_PI)*(cos(M_PI*y1/coreLength) - 
			cos(M_PI*y2/coreLength));
			double sx = (1./model.dx)*(xLength/M_PI)*(cos(M_PI*x1/xLength) - 
			cos(M_PI*x2/xLength));

			if (y < coreLength){
				s = 1.e13*sx*sy;
			}
			else{
				s = 0.0;
			}
			model.setCellNeutronFlux(i, j, s);
		}
	}

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);

	// Sets the flux limiter type
	spec.setFluxLimiter(limiter);

	// Adds the speices
	ids = spec.addSpeciesFromFile(speciesNamesFile);
	std::vector<double> bcs = {};

	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	readCSV(refSolData, std::string(path + "solutionMSR2DDepletionSmall.csv"));

	// Print species name to file
	if (myid == 0){
		pOutputFile = fopen(outputFileName.c_str(), "a");
		fprintf(pOutputFile, "%s %s %s ", "variables", "x", "y");
		for (int id; id < ids.size(); id++){
			std::string name = spec.getSpeciesName(id);
			fprintf(pOutputFile, "%s ", name.c_str());
		}
		fprintf(pOutputFile, "\n");
		fclose(pOutputFile);
	}

	// Sets the species sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);

	// Adds the boundary conditions
	spec.setBoundaryCondition("newmann","east", ids, bcs);
	spec.setBoundaryCondition("newmann","west", ids, bcs);
	spec.setBoundaryCondition("periodic","south", ids);
	spec.setBoundaryCondition("periodic","north", ids);

	// Writes transition matrix and initial condition
	if (myid == 0){
		spec.writeTransitionMatrixToFile("transitionMatrixMSR2DDepletion.csv");
		spec.writeInitialConditionToFile("initialConditionMSR2DDepletion.csv");
		pOutputFile = fopen(outputFileName.c_str(), "a");
	}

	// Loops over the time steps to solve
	double totalSolveTime = 0;
	for (int k = 0; k < steps; k++){
		auto start = std::chrono::high_resolution_clock::now();
		t = t + dt;
		//std::cout << "solving" << std::endl;
		spec.solve(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		if (myid == 0){
			solveTime = duration.count()/1.e6;
			totalSolveTime += solveTime;

			fprintf(pOutputFile, "time: %5.4e\n", t);
			int index = 0;
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					fprintf(pOutputFile, "%5.4e %5.4e ", xc, yc);

					for (int id = 0; id < ids.size(); id++){
						std::string name = spec.getSpeciesName(ids[id]);
						double con = spec.getSpecies(i, j, ids[id]);
						fprintf(pOutputFile, "%17.16e ", con);
						solData(index, k) = con;
						index ++;
					}
					fprintf(pOutputFile, "\n");
				}
			}
			VectorD refSolVect = refSolData.col(k), solVect = solData.col(k);
			//double rmse = computeRelativeRMSE(refSolVect, solVect);
			double Einf = computeRelativeEinfty(refSolVect, solVect);
			double rmse = 0;
			printf("Solve Step: %d %s Solve Time: %f Einf %e\n", k, solverType.c_str(), solveTime, Einf);
		}
	}
	if (myid == 0){
		double Einf = computeRelativeEinfty(refSolData, solData);
		printf("%s Solve Time: %f Einf %e\n", solverType.c_str(), totalSolveTime, Einf);
		fprintf(pOutputFile, "end\n");
		//writeCSV(solData, outputFileNameMatrix);
	}
}

//*****************************************************************************
// MSR lump depletion
//
// small nuclide list
//*****************************************************************************
void msrDepletionSmallLumped(int myid, std::string solverType){
	double t = 0.0;
	//int steps = 1;
	int steps = 10;
	//double depletionTime = 20.; // Days
	double depletionTime = 200.; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 1, yCells = 1;
	double xLength = 1.92024, yLength = 0.6858, xc, yc;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string solPath = getDataPath() + "caseStudy/msrLumpDepletion/";
	std::string outputFileName = solverType+"Steps0MSRDepletionSmallLumped.out";
	std::string outputFileNameMatrix = solverType+"MSRDepletionSmallLumped.csv";
	std::string limiter = "First order upwind";
	std::ofstream outputFile;
	outputFile.open(outputFileName, std::ios::trunc);
	MatrixD refSolData;
	FILE * pOutputFile;
	meshCell* cell = nullptr;
	pOutputFile = fopen(outputFileName.c_str(), "a");
	std::vector<std::string> gasNames = {"H-1Gas", "N-14Gas", "O-16Gas", "Kr-83Gas", 
													 "Xe-131Gas", "Xe-135Gas"};

	//printf("Time step size: %f\n", totalTime/(double)steps);

	if (myid == 0){
		fprintf(pOutputFile, "Solver: %s \n", solverType.c_str());
		fprintf(pOutputFile, "Total problem time: %4.3e\n", totalTime);
		fprintf(pOutputFile, "yLength: %4.3e\n", yLength);
		fprintf(pOutputFile, "xLength: %4.3e\n", xLength);
		fprintf(pOutputFile, "dt: %8.7e\n", dt);
		fprintf(pOutputFile, "yCells: %d\n", yCells);
		fprintf(pOutputFile, "xCells: %d\n", xCells);
		fprintf(pOutputFile, "%s %s %s %s %s %s %s %s %s \n", "variables", 
			"i", "j", "x", "y", "Name", "Con", "solCon", "Einf");
	fclose(pOutputFile);
	}

	// Sets file paths for the input data	
	std::string speciesNamesFile = path + "speciesInputNamesMSRMassTransportSmall.dat";
	std::string speciesDecayFile = path + "speciesInputDecayMSRMassTransportSmall.dat";
	std::string speciesTransFile = path + "speciesInputTransMSRMassTransportSmall.dat";

	// Builds the model
	modelMesh model(xCells, yCells, xLength, yLength);

	// Builds the species object
	speciesDriver spec = speciesDriver(&model);

	// Sets the neutron flux
	model.setSystemNeutronFlux(4.05283329E+12);
	model.setSystemGasInterfacialAreaCon(2.42994803E+01);
	model.setSystemTemperature(912.);
	model.setSystemGasVoidFraction(3.00E-04);
	model.setSystemWallInterfacialAreaCon(1.35925218E+02);

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);
	//spec.setIntegratorSolver("implicit", solverType);

	// Adds the speices
	ids = spec.addSpeciesFromFile(speciesNamesFile);

	// Sets the species sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);
	spec.setGasSpargingFromFile(path + "speciesInputGasSpargingSmall.txt");
	spec.setWallDepositionFromFile(path + "speciesInputWallDepositionSmall.txt");

	// Sets the gas removal
	for (int k = 0; k < gasNames.size(); k++){
		spec.setRemoval(0, 0, spec.getSpeciesID(gasNames[k]), -2.05e-4);
	}

	// Writes transition matrix and initial condition
	//spec.writeTransitionMatrixToFile("transitionMatrixMSRLumpDepletionSmall.csv");
	//spec.writeInitialConditionToFile("initialConditionMSRLumpDepletionSmall.csv");
	pOutputFile = fopen(outputFileName.c_str(), "a");

	//// Gets the solution
	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	readCSV(refSolData, std::string(solPath + "solutionMSRLumpDepletionSmall.csv"));

	// Loops over the time steps to solve
	double totalSolveTime = 0;
	for (int k = 0; k < steps; k++){
		auto start = std::chrono::high_resolution_clock::now();
		t = t + dt;
		spec.solve(t);
		//spec.solveImplicit(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		solveTime = duration.count()/1.e6;
		totalSolveTime += solveTime;

		if (myid == 0){

			cell = model.getCellByLoc(0,0); xc = cell->x; yc = cell->y;

			fprintf(pOutputFile, "time: %5.4e\n", t);
			//printf("time: %5.4e\n", t);

			for (int id = 0; id < ids.size(); id++){
				std::string name = spec.getSpeciesName(ids[id]);
				double con = spec.getSpecies(0, 0, ids[id]);
				double solCon = refSolData(id,k);
				//fprintf(pOutputFile, "%5.4e %5.4e %s %17.16e \n", xc, yc, name.c_str(), con);
				fprintf(pOutputFile, "%d %d %5.4e %5.4e %s %17.16e %17.16e %8.7e\n", 0, 0, xc, yc, name.c_str(), con, solCon, abs(con-solCon)/solCon);
				//if (k == steps -1){
				//	printf("%s %17.16e %17.16e %8.7e\n", name.c_str(), con, solCon, abs(con-solCon)/solCon);
				//}
				solData(id, k) = con;
			}
			fprintf(pOutputFile, "\n");
			VectorD refSol = refSolData.col(k), sol = solData.col(k);
			//VectorD refSol = refSolData.col(0), sol = solData.col(k);
			double E1 = computeRelativeE1(refSol, sol);
			double E2 = computeRelativeE2(refSol, sol);
			double Einf = computeRelativeEinfty(refSol, sol);
			//double E1 = 0.;
			//double E2 = 0.;
			//double Einf  = 0.;
			printf("%s %e %e %e %e\n", solverType.c_str(), t, Einf, E1, E2);
			//printf("%e %e %e\n", Einf, E1, E2);
		}
	}
	if (myid == 0){
		//VectorD refSol = refSolData.col(0), sol = solData.col(steps-1);
		//double E1 = computeRelativeE1(refSol, sol);
		//double E2 = computeRelativeE2(refSol, sol);
		//double Einf = computeRelativeEinfty(refSol, sol);
		fprintf(pOutputFile, "end\n");
		printf("%s Solve Time: %f \n", solverType.c_str(), totalSolveTime);
		//printf("%s Solve Time: %f %e %e %e \n", solverType.c_str(), 
		//	totalSolveTime, Einf, E1, E2);
		//writeCSV(solData, outputFileNameMatrix);
	}
}

//*****************************************************************************
// MSR lump depletion
//
// medium nuclide list
//*****************************************************************************
void msrDepletionMediumLumped(int myid, std::string solverType){
	double t = 0.0;
	//int steps = 100000;
	int steps = 10;
	double depletionTime = 200.; // Days
	//double depletionTime = 20.; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 1, yCells = 1;
	double xLength = 1.92024, yLength = 0.6858, xc, yc;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string solPath = getDataPath() + "caseStudy/msrLumpDepletion/";
	std::string outputFileName = solverType+"Steps0MSRDepletionMediumLumped.out";
	std::string outputFileNameMatrix = solverType+"MSRDepletionMediumLumped.csv";
	std::string limiter = "First order upwind";
	std::ofstream outputFile;
	outputFile.open(outputFileName, std::ios::trunc);
	MatrixD refSolData;
	FILE * pOutputFile;
	meshCell* cell = nullptr;
	pOutputFile = fopen(outputFileName.c_str(), "a");
	std::vector<std::string> gasNames = {"H-1Gas", "N-14Gas", "O-16Gas", "Kr-83Gas", 
													 "Xe-131Gas", "Xe-133Gas", "Xe-135Gas"};

	//printf("Time step size: %f\n", totalTime/(double)steps);

	if (myid == 0){
		fprintf(pOutputFile, "Solver: %s \n", solverType.c_str());
		fprintf(pOutputFile, "Total problem time: %4.3e\n", totalTime);
		fprintf(pOutputFile, "yLength: %4.3e\n", yLength);
		fprintf(pOutputFile, "xLength: %4.3e\n", xLength);
		fprintf(pOutputFile, "dt: %8.7e\n", dt);
		fprintf(pOutputFile, "yCells: %d\n", yCells);
		fprintf(pOutputFile, "xCells: %d\n", xCells);
		fprintf(pOutputFile, "%s %s %s %s %s %s %s %s %s \n", "variables", 
			"i", "j", "x", "y", "Name", "Con", "solCon", "Einf");
		fclose(pOutputFile);
	}

	// Sets file paths for the input data	
	std::string speciesNamesFile = path + "speciesInputNamesMSRMassTransportMedium.dat";
	std::string speciesDecayFile = path + "speciesInputDecayMSRMassTransportMedium.dat";
	std::string speciesTransFile = path + "speciesInputTransMSRMassTransportMedium.dat";

	// Builds the model
	modelMesh model(xCells, yCells, xLength, yLength);

	// Builds the species object
	speciesDriver spec = speciesDriver(&model);

	// Sets the neutron flux
	model.setSystemNeutronFlux(4.05283329E+12);
	model.setSystemGasInterfacialAreaCon(2.42994803E+01);
	model.setSystemTemperature(912.);
	model.setSystemGasVoidFraction(3.00E-04);
	model.setSystemWallInterfacialAreaCon(1.35925218E+02);

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);
	//spec.setIntegratorSolver("implicit", solverType);

	// Adds the speices
	ids = spec.addSpeciesFromFile(speciesNamesFile);

	// Sets the species sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);
	spec.setGasSpargingFromFile(path + "speciesInputGasSpargingMedium.txt");
	spec.setWallDepositionFromFile(path + "speciesInputWallDepositionMedium.txt");

	// Sets the gas removal
	for (int k = 0; k < gasNames.size(); k++){
		spec.setRemoval(0, 0, spec.getSpeciesID(gasNames[k]), -2.05e-4);
	}

	// Writes transition matrix and initial condition
	spec.writeTransitionMatrixToFile("transitionMatrixMSRLumpDepletionMedium.csv");
	spec.writeInitialConditionToFile("initialConditionMSRLumpDepletionMedium.csv");
	pOutputFile = fopen(outputFileName.c_str(), "a");

	// Gets the solution
	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	readCSV(refSolData, std::string(solPath + "solutionMSRLumpDepletionMedium.csv"));

	// Loops over the time steps to solve
	double totalSolveTime = 0;
	for (int k = 0; k < steps; k++){
		auto start = std::chrono::high_resolution_clock::now();
		t = t + dt;
		spec.solve(t);
		//spec.solveImplicit(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		solveTime = duration.count()/1.e6;
		totalSolveTime += solveTime;
		if (myid == 0){

			cell = model.getCellByLoc(0,0); xc = cell->x; yc = cell->y;

			fprintf(pOutputFile, "time: %5.4e\n", t);
			//printf("time: %5.4e\n", t);

			for (int id = 0; id < ids.size(); id++){
				std::string name = spec.getSpeciesName(ids[id]);
				double con = spec.getSpecies(0, 0, ids[id]);
				double solCon = refSolData(id,k);
				fprintf(pOutputFile, "%d %d %5.4e %5.4e %s %17.16e %17.16e %8.7e\n", 0, 0, xc, yc, name.c_str(), con, solCon, abs(con-solCon)/solCon);
				//printf("%5.4e %5.4e %s %17.16e \n", xc, yc, name.c_str(), con);
				solData(id, k) = con;
			}
			fprintf(pOutputFile, "\n");
			VectorD refSol = refSolData.col(k), sol = solData.col(k);
			double E1 = computeRelativeE1(refSol, sol);
			double E2 = computeRelativeE2(refSol, sol);
			double Einf = computeRelativeEinfty(refSol, sol);
			//double E1 = 0.;
			//double E2 = 0.;
			//double Einf  = 0.;
			printf("%e %e %e\n", Einf, E1, E2);
		}
	}
	if (myid == 0){
		//VectorD refSol = refSolData.col(0), sol = solData.col(steps-1);
		//double E1 = computeRelativeE1(refSol, sol);
		//double E2 = computeRelativeE2(refSol, sol);
		//double Einf = computeRelativeEinfty(refSol, sol);
		fprintf(pOutputFile, "end\n");
		printf("%s Solve Time: %f \n", solverType.c_str(), totalSolveTime);
		//printf("%s Solve Time: %f %e %e %e \n", solverType.c_str(), 
		//	totalSolveTime, Einf, E1, E2);
		//writeCSV(solData, outputFileNameMatrix);
	}
}
//*****************************************************************************
// MSR 2-D depletion
//
// Small nuclide list case study for the 3x9 case
//*****************************************************************************
void msr2DDepletionSmall3x9(int myid, std::string solverType){
	double t = 0.0;
	int steps = 10;
	double coreLength = 1.92024;
	double xLength = 0.6858, yLength = 3.*coreLength;		// Meters
	double v_y = 0.25;											// m/s
	double depletionTime = 200.;							// Days
	//double depletionTime = 20.;							// Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 3, yCells = 9;
	double xc, yc, s;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string solPath = getDataPath() + "caseStudy/2DMSRTransport/";
	std::string outputFileName = solverType+"Steps0MSRDepletionSmall3x9.out";
	std::string outputFileNameMatrix = solverType+"MSR2DDepletionSmall3x9.csv";
	std::string limiter = "First order upwind";
	std::ofstream outputFile;
	outputFile.open(outputFileName, std::ios::trunc);
	MatrixD refSolData;
	FILE * pOutputFile;
	meshCell* cell = nullptr;
	std::vector<std::string> wallTransNames={"Nb-93", "Mo-95", "Tc-99", "Ru-106", "Rh-103",
													     "Rh-105", "Ag-109", "Sn-126"};
	std::vector<std::string> gasNames = {"H-1Gas", "N-14Gas", "O-16Gas", "Kr-83Gas", 
													 "Xe-131Gas", "Xe-135Gas"};
	// {graphite, heat exchanger, piping}
	std::vector<double> massTransCoeffs = {5.334e-6, 4.657e-5, 1.041e-4};

	if (myid == 0){
		pOutputFile = fopen(outputFileName.c_str(), "a");
		fprintf(pOutputFile, "Solver: %s \n", solverType.c_str());
		fprintf(pOutputFile, "Total problem time: %4.3e\n", totalTime);
		fprintf(pOutputFile, "yLength: %4.3e\n", yLength);
		fprintf(pOutputFile, "xLength: %4.3e\n", xLength);
		fprintf(pOutputFile, "dt: %8.7e\n", dt);
		fprintf(pOutputFile, "yCells: %d\n", yCells);
		fprintf(pOutputFile, "xCells: %d\n", xCells);
		fprintf(pOutputFile, "%s %s %s %s %s %s %s %s %s \n", "variables", 
			"i", "j", "x", "y", "Name", "Con", "solCon", "Einf");
		fclose(pOutputFile);
	}

	// Sets file paths for the input data	
	std::string speciesNamesFile = path + "speciesInputNamesMSRMassTransportSmall.dat";
	std::string speciesDecayFile = path + "speciesInputDecayMSRMassTransportSmall.dat";
	std::string speciesTransFile = path + "speciesInputTransMSRMassTransportSmall.dat";

	// Builds the model
	modelMesh model(xCells, yCells, xLength, yLength);

	// Adds boundary surface
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");

	// set y velocity
	model.setConstantYVelocity(v_y);

	// Sets the variables
	model.setSystemTemperature(path + "temperatureMap3x9.csv");
	model.setSystemNeutronFlux(path + "neutronFluxMap3x9.csv");
	model.setSystemGasVoidFraction(path + "voidMap3x9.csv");
	model.setSystemGasInterfacialAreaCon(path + "gasIntAreaMap3x9.csv");
	model.setSystemWallInterfacialAreaCon(path + "wallIntAreaMap3x9.csv");

	// Builds the species object
	speciesDriver spec = speciesDriver(&model);

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);
	//spec.setIntegratorSolver("implicit", solverType);

	// Sets the flux limiter type
	spec.setFluxLimiter(limiter);

	// Adds the speices
	ids = spec.addSpeciesFromFile(speciesNamesFile);
	spec.setGasSpargingFromFile(path + "speciesInputGasSpargingSmall.txt");

	// Loops through to build the liqIDs and wallIDs for wall deposition 
	std::vector<int> liqIDs = {}, wallIDs = {};
	std::vector<bool> infSinks = {};
	for (int k = 0; k < wallTransNames.size(); k++){
		std::string name = wallTransNames[k];
		liqIDs.push_back(spec.getSpeciesID(name+"Liq"));
		wallIDs.push_back(spec.getSpeciesID(name+"Wall"));
		infSinks.push_back(true);
	}

	// Loops through the model to set the wall mass transport coefficients
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			// The core region. This contains graphite
			if (j > 0 and j < 3){
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[0]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}
			// Heat exchanger
			else if (j == 4){
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[1]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}
			else{
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[2]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}

			// Adds the gas removal
			if (j == 3){
				for (int k = 0; k < gasNames.size(); k++){
					spec.setRemoval(i, j, spec.getSpeciesID(gasNames[k]), -2.05e-4);
				}
			}

		}
	}

	// Boundary condition vector
	std::vector<double> bcs = {};

	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	readCSV(refSolData, std::string(solPath + "solutionMSR2DDepletionSmall3x9.csv"));

	// Print species name to file
	//if (myid == 0){
	//	pOutputFile = fopen(outputFileName.c_str(), "a");
	//	fprintf(pOutputFile, "%s %s %s ", "variables", "x", "y");
	//	for (int id; id < ids.size(); id++){
	//		std::string name = spec.getSpeciesName(id);
	//		fprintf(pOutputFile, "%s ", name.c_str());
	//	}
	//	fprintf(pOutputFile, "\n");
	//	fclose(pOutputFile);
	//}

	// Sets the species sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);

	// Adds the boundary conditions
	spec.setBoundaryCondition("newmann","east", ids, bcs);
	spec.setBoundaryCondition("newmann","west", ids, bcs);
	spec.setBoundaryCondition("periodic","south", ids);
	spec.setBoundaryCondition("periodic","north", ids);

	// Writes transition matrix and initial condition
	if (myid == 0){
		spec.writeTransitionMatrixToFile("transitionMatrixMSR2DDepletionSmall3x9.csv");
		spec.writeInitialConditionToFile("initialConditionMSR2DDepletionSmall3x9.csv");
	}

	// Loops over the time steps to solve
	double totalSolveTime = 0;
	for (int k = 0; k < steps; k++){
		auto start = std::chrono::high_resolution_clock::now();
		t = t + dt;
		spec.solve(t);
		//spec.solveImplicit(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		if (myid == 0){
			pOutputFile = fopen(outputFileName.c_str(), "a");
			solveTime = duration.count()/1.e6;
			totalSolveTime += solveTime;

			fprintf(pOutputFile, "time: %5.4e\n", t);
			int index = 0;
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j); xc = cell->x; yc = cell->y;

					for (int id = 0; id < ids.size(); id++){
						std::string name = spec.getSpeciesName(ids[id]);
						double con = spec.getSpecies(i, j, ids[id]);
						double solCon = refSolData(index,k);
						fprintf(pOutputFile, "%d %d %5.4e %5.4e %s %17.16e %17.16e %8.7e\n", i, j, xc, yc, name.c_str(), con, solCon, abs(con-solCon)/solCon);
						solData(index, k) = con;
						index ++;
					}
				}
			}
			fprintf(pOutputFile, "\n");
			VectorD refSolVect = refSolData.col(k), solVect = solData.col(k);
			double E1 = computeRelativeE1(refSolVect, solVect);
			double E2 = computeRelativeE2(refSolVect, solVect);
			double Einf = computeRelativeEinfty(refSolVect, solVect);
			//double E1 = 0.;
			//double E2 = 0.;
			//double Einf  = 0.;
			printf("%e %e %e\n", Einf, E1, E2);
			fclose(pOutputFile);
		}
	}
	if (myid == 0){
		pOutputFile = fopen(outputFileName.c_str(), "a");
		//VectorD refSol = refSolData.col(0), sol = solData.col(steps-1);
		//double E1 = computeRelativeE1(refSol, sol);
		//double E2 = computeRelativeE2(refSol, sol);
		//double Einf = computeRelativeEinfty(refSol, sol);
		printf("%s Solve Time: %f \n", solverType.c_str(), totalSolveTime);
		//printf("%s Solve Time: %f %e %e %e \n", solverType.c_str(), 
		//	totalSolveTime, Einf, E1, E2);
		fprintf(pOutputFile, "end\n");
		fclose(pOutputFile);
		//writeCSV(solData, outputFileNameMatrix);
	}
}

//*****************************************************************************
// MSR 2-D depletion
//
// Medium nuclide list case study for the 3x9 case
//*****************************************************************************
void msr2DDepletionMedium3x9(int myid, std::string solverType){
	double t = 0.0;
	int steps = 10;
	double coreLength = 1.92024;
	double xLength = 0.6858, yLength = 3.*coreLength;		// Meters
	double v_y = 0.25;											// m/s
	double depletionTime = 200.;							// Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 3, yCells = 9;
	double xc, yc, s;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string solPath = getDataPath() + "caseStudy/2DMSRTransport/";
	std::string outputFileName = "msr2DDepletionMedium3x9.out";
	std::string outputFileNameMatrix = solverType+"MSR2DDepletionMedium3x9.csv";
	std::string limiter = "First order upwind";
	std::ofstream outputFile;
	outputFile.open(outputFileName, std::ios_base::app);
	MatrixD refSolData;
	FILE * pOutputFile;
	meshCell* cell = nullptr;
	std::vector<std::string> wallTransNames={"Nb-93", "Mo-95", "Mo-97", "Mo-98", "Mo-99",
													     "Mo-100", "Tc-99", "Ru-101", "Ru-102", "Ru-103", 
														  "Ru-104", "Ru-106", "Rh-105", "Rh-103", "Pd-105",
														  "Pd-107", "Pd-108", "Ag-109", "Sn-126"};
	std::vector<std::string> gasNames = {"H-1Gas", "N-14Gas", "O-16Gas", "Kr-83Gas", 
													 "Xe-131Gas", "Xe-133Gas", "Xe-135Gas"};
	// {graphite, heat exchanger, piping}
	std::vector<double> massTransCoeffs = {5.334e-6, 4.657e-5, 1.041e-4};

	if (myid == 0){
		pOutputFile = fopen(outputFileName.c_str(), "a");
		fprintf(pOutputFile, "Solver: %s \n", solverType.c_str());
		fprintf(pOutputFile, "Total problem time: %4.3e\n", totalTime);
		fprintf(pOutputFile, "yLength: %4.3e\n", yLength);
		fprintf(pOutputFile, "xLength: %4.3e\n", xLength);
		fprintf(pOutputFile, "dt: %8.7e\n", dt);
		fprintf(pOutputFile, "yCells: %d\n", yCells);
		fprintf(pOutputFile, "xCells: %d\n", xCells);
		fclose(pOutputFile);
	}

	// Sets file paths for the input data	
	std::string speciesNamesFile = path + "speciesInputNamesMSRMassTransportMedium.dat";
	std::string speciesDecayFile = path + "speciesInputDecayMSRMassTransportMedium.dat";
	std::string speciesTransFile = path + "speciesInputTransMSRMassTransportMedium.dat";

	// Builds the model
	modelMesh model(xCells, yCells, xLength, yLength);

	// Adds boundary surface
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");

	// set y velocity
	model.setConstantYVelocity(v_y);

	// Builds the species object
	speciesDriver spec = speciesDriver(&model);

	// Sets the variables
	model.setSystemTemperature(path + "temperatureMap3x9.csv");
	model.setSystemNeutronFlux(path + "neutronFluxMap3x9.csv");
	model.setSystemGasVoidFraction(path + "voidMap3x9.csv");
	model.setSystemGasInterfacialAreaCon(path + "gasIntAreaMap3x9.csv");
	model.setSystemWallInterfacialAreaCon(path + "wallIntAreaMap3x9.csv");

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);

	// Sets the flux limiter type
	spec.setFluxLimiter(limiter);

	// Adds the speices
	ids = spec.addSpeciesFromFile(speciesNamesFile);
	spec.setGasSpargingFromFile(path + "speciesInputGasSpargingMedium.txt");

	// Loops through to build the liqIDs and wallIDs for wall deposition 
	std::vector<int> liqIDs = {}, wallIDs = {};
	std::vector<bool> infSinks = {};
	for (int k = 0; k < wallTransNames.size(); k++){
		std::string name = wallTransNames[k];
		liqIDs.push_back(spec.getSpeciesID(name+"Liq"));
		wallIDs.push_back(spec.getSpeciesID(name+"Wall"));
		infSinks.push_back(true);
	}

	// Loops through the model to set the wall mass transport coefficients
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			// The core region. This contains graphite
			if (j > 0 and j < 3){
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[0]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}
			// Heat exchanger
			else if (j == 4){
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[1]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}
			else{
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[2]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}

			// Adds the gas removal
			if (j == 3){
				for (int k = 0; k < gasNames.size(); k++){
					spec.setRemoval(i, j, spec.getSpeciesID(gasNames[k]), -2.05e-4);
				}
			}
		}
	}

	// boundary condition vector
	std::vector<double> bcs = {};

	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	readCSV(refSolData, std::string(solPath + "solutionMSR2DDepletionMedium3x9.csv"));

	// Print species name to file
	if (myid == 0){
		pOutputFile = fopen(outputFileName.c_str(), "a");
		fprintf(pOutputFile, "%s %s %s ", "variables", "x", "y");
		for (int id; id < ids.size(); id++){
			std::string name = spec.getSpeciesName(id);
			fprintf(pOutputFile, "%s ", name.c_str());
		}
		fprintf(pOutputFile, "\n");
		fclose(pOutputFile);
	}

	// Sets the species sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);

	// Adds the boundary conditions
	spec.setBoundaryCondition("newmann","east", ids, bcs);
	spec.setBoundaryCondition("newmann","west", ids, bcs);
	spec.setBoundaryCondition("periodic","south", ids);
	spec.setBoundaryCondition("periodic","north", ids);

	// Writes transition matrix and initial condition
	if (myid == 0){
		spec.writeTransitionMatrixToFile("transitionMatrixMSR2DDepletionMedium3x9.csv");
		spec.writeInitialConditionToFile("initialConditionMSR2DDepletionMedium3x9.csv");
		pOutputFile = fopen(outputFileName.c_str(), "a");
	}

	// Loops over the time steps to solve
	double totalSolveTime = 0;
	for (int k = 0; k < steps; k++){
		auto start = std::chrono::high_resolution_clock::now();
		t = t + dt;
		//std::cout << "solving" << std::endl;
		spec.solve(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		if (myid == 0){
			solveTime = duration.count()/1.e6;
			totalSolveTime += solveTime;

			fprintf(pOutputFile, "time: %5.4e\n", t);
			int index = 0;
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					cell = model.getCellByLoc(i,j); xc = cell->x; yc = cell->y;
					fprintf(pOutputFile, "%5.4e %5.4e ", xc, yc);
					cell = model.getCellByLoc(i,j); xc = cell->x; yc = cell->y;
					double temp = cell->getTemperature(), fract = cell->getGasVoidFraction();
					double flux = cell->getNeutronFlux(), wallInt = cell->getWallInterfacialAreaCon();
					double gasInt = cell->getGasInterfacialAreaCon();
					//printf("%i %i %5.4e %5.4e %8.7e %8.7e %8.7e %8.7e %8.7e \n", 
					//	i, j, xc, yc, temp, fract, flux, wallInt, gasInt);

					for (int id = 0; id < ids.size(); id++){
						std::string name = spec.getSpeciesName(ids[id]);
						double con = spec.getSpecies(i, j, ids[id]);
						fprintf(pOutputFile, "%17.16e ", con);
						//printf("%5.4e %5.4e %s %17.16e \n", xc, yc, name.c_str(), con);
						solData(index, k) = con;
						index ++;
					}
					fprintf(pOutputFile, "\n");
				}
			}
			VectorD refSolVect = refSolData.col(k), solVect = solData.col(k);
			double E1 = computeRelativeE1(refSolVect, solVect);
			double E2 = computeRelativeE2(refSolVect, solVect);
			double Einf = computeRelativeEinfty(refSolVect, solVect);
			//double E1 = 0.;
			//double E2 = 0.;
			//double Einf  = 0.;
			printf("%s %e %e %e %e\n", solverType.c_str(), t, Einf, E1, E2);
		}
	}
	if (myid == 0){
		printf("%s Solve Time: %f\n", solverType.c_str(), totalSolveTime);
		fprintf(pOutputFile, "end\n");
		writeCSV(solData, outputFileNameMatrix);
	}
}

//*****************************************************************************
// MSR 2-D depletion
//
// Small nuclide list case study for the 9x27 case
//*****************************************************************************
void msr2DDepletionSmall9x27(int myid, std::string solverType){
	double t = 0.0;
	int steps = 10;
	double coreLength = 1.92024;
	double xLength = 0.6858, yLength = 3.*coreLength;		// Meters
	double v_y = 0.25;											// m/s
	double depletionTime = 200.;							// Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 9, yCells = 27;
	double xc, yc, s;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string solPath = getDataPath() + "caseStudy/msrLumpDepletion/";
	std::string outputFileName = "msr2DDepletionSmall9x27.out";
	std::string outputFileNameMatrix = solverType+"MSR2DDepletionSmall9x27.csv";
	std::string limiter = "First order upwind";
	std::ofstream outputFile;
	outputFile.open(outputFileName, std::ios_base::app);
	MatrixD refSolData;
	FILE * pOutputFile;
	meshCell* cell = nullptr;
	std::vector<std::string> wallTransNames={"Nb-93", "Mo-95", "Tc-99", "Ru-106", "Rh-103",
													     "Rh-105", "Ag-109", "Sn-126"};
	std::vector<std::string> gasNames = {"H-1Gas", "N-14Gas", "O-16Gas", "Kr-83Gas", 
													 "Xe-131Gas", "Xe-135Gas"};
	// {graphite, heat exchanger, piping}
	std::vector<double> massTransCoeffs = {5.334e-6, 4.657e-5, 1.041e-4};

	if (myid == 0){
		pOutputFile = fopen(outputFileName.c_str(), "a");
		fprintf(pOutputFile, "Solver: %s \n", solverType.c_str());
		fprintf(pOutputFile, "Total problem time: %4.3e\n", totalTime);
		fprintf(pOutputFile, "yLength: %4.3e\n", yLength);
		fprintf(pOutputFile, "xLength: %4.3e\n", xLength);
		fprintf(pOutputFile, "dt: %8.7e\n", dt);
		fprintf(pOutputFile, "yCells: %d\n", yCells);
		fprintf(pOutputFile, "xCells: %d\n", xCells);
		fclose(pOutputFile);
	}

	// Sets file paths for the input data	
	std::string speciesNamesFile = path + "speciesInputNamesMSRMassTransportSmall.dat";
	std::string speciesDecayFile = path + "speciesInputDecayMSRMassTransportSmall.dat";
	std::string speciesTransFile = path + "speciesInputTransMSRMassTransportSmall.dat";

	// Builds the model
	modelMesh model(xCells, yCells, xLength, yLength);

	// Adds boundary surface
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");

	// set y velocity
	model.setConstantYVelocity(v_y);

	// Builds the species object
	speciesDriver spec = speciesDriver(&model);

	// Sets the variables
	model.setSystemTemperature(path + "temperatureMap9x27.csv");
	model.setSystemNeutronFlux(path + "neutronFluxMap9x27.csv");
	model.setSystemGasVoidFraction(path + "voidMap9x27.csv");
	model.setSystemGasInterfacialAreaCon(path + "gasIntAreaMap9x27.csv");
	model.setSystemWallInterfacialAreaCon(path + "wallIntAreaMap9x27.csv");

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);

	// Sets the flux limiter type
	spec.setFluxLimiter(limiter);

	// Adds the speices
	ids = spec.addSpeciesFromFile(speciesNamesFile);
	spec.setGasSpargingFromFile(path + "speciesInputGasSpargingSmall.txt");

	// Loops through to build the liqIDs and wallIDs for wall deposition 
	std::vector<int> liqIDs = {}, wallIDs = {};
	std::vector<bool> infSinks = {};
	for (int k = 0; k < wallTransNames.size(); k++){
		std::string name = wallTransNames[k];
		liqIDs.push_back(spec.getSpeciesID(name+"Liq"));
		wallIDs.push_back(spec.getSpeciesID(name+"Wall"));
		infSinks.push_back(true);
	}

	// Loops through the model to set the wall mass transport coefficients
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			// The core region. This contains graphite
			if (j > 0 and j < 9){
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[0]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}
			// Heat exchanger
			else if (j >= 12 and j <= 14){
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[1]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}
			else{
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[2]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}

			// Adds the gas removal
			if (j == 10){
				for (int k = 0; k < gasNames.size(); k++){
					spec.setRemoval(i, j, spec.getSpeciesID(gasNames[k]), -2.05e-4);
				}
			}
		}
	}

	// Sets the boundary condition values
	std::vector<double> bcs = {};

	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	//readCSV(refSolData, std::string(path + "solutionMSR2DDepletionSmall.csv"));

	// Print species name to file
	if (myid == 0){
		pOutputFile = fopen(outputFileName.c_str(), "a");
		fprintf(pOutputFile, "%s %s %s ", "variables", "x", "y");
		for (int id; id < ids.size(); id++){
			std::string name = spec.getSpeciesName(id);
			fprintf(pOutputFile, "%s ", name.c_str());
		}
		fprintf(pOutputFile, "\n");
		fclose(pOutputFile);
	}

	// Sets the species sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);

	// Adds the boundary conditions
	spec.setBoundaryCondition("newmann","east", ids, bcs);
	spec.setBoundaryCondition("newmann","west", ids, bcs);
	spec.setBoundaryCondition("periodic","south", ids);
	spec.setBoundaryCondition("periodic","north", ids);

	// Writes transition matrix and initial condition
	if (myid == 0){
		spec.writeTransitionMatrixToFile("transitionMatrixMSR2DDepletionSmall9x27.csv");
		spec.writeInitialConditionToFile("initialConditionMSR2DDepletionSmall9x27.csv");
		pOutputFile = fopen(outputFileName.c_str(), "a");
	}

	// Loops over the time steps to solve
	double totalSolveTime = 0;
	for (int k = 0; k < steps; k++){
		auto start = std::chrono::high_resolution_clock::now();
		t = t + dt;
		spec.solve(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		if (myid == 0){
			solveTime = duration.count()/1.e6;
			totalSolveTime += solveTime;

			fprintf(pOutputFile, "time: %5.4e\n", t);
			int index = 0;
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					fprintf(pOutputFile, "%5.4e %5.4e ", xc, yc);
					cell = model.getCellByLoc(i,j); xc = cell->x; yc = cell->y;
					double temp = cell->getTemperature(), fract = cell->getGasVoidFraction();
					double flux = cell->getNeutronFlux(), wallInt = cell->getWallInterfacialAreaCon();
					double gasInt = cell->getGasInterfacialAreaCon();
					printf("%i %i %5.4e %5.4e %8.7e %8.7e %8.7e %8.7e %8.7e \n", 
						i, j, xc, yc, temp, fract, flux, wallInt, gasInt);

					for (int id = 0; id < ids.size(); id++){
						std::string name = spec.getSpeciesName(ids[id]);
						double con = spec.getSpecies(i, j, ids[id]);
						fprintf(pOutputFile, "%17.16e ", con);
						printf("%5.4e %5.4e %s %17.16e \n", xc, yc, name.c_str(), con);
						solData(index, k) = con;
						index ++;
					}
					fprintf(pOutputFile, "\n");
				}
			}
			//VectorD refSolVect = refSolData.col(k), solVect = solData.col(k);
			//double E1 = computeRelativeE1(refSol, sol);
			//double E2 = computeRelativeE2(refSol, sol);
			//double Einf = computeRelativeEinfty(refSol, sol);
			double E1 = 0.;
			double E2 = 0.;
			double Einf  = 0.;
			printf("%s %e %e %e %e\n", solverType.c_str(), t, Einf, E1, E2);
		}
	}
	if (myid == 0){
		printf("%s Solve Time: %f\n", solverType.c_str(), totalSolveTime);
		fprintf(pOutputFile, "end\n");
		writeCSV(solData, outputFileNameMatrix);
	}
}

//*****************************************************************************
// MSR 2-D depletion
//
// Medium nuclide list case study for the 9x27 case
//*****************************************************************************
void msr2DDepletionMedium9x27(int myid, std::string solverType){
	double t = 0.0;
	int steps = 10;
	double coreLength = 1.92024;
	double xLength = 0.6858, yLength = 3.*coreLength;		// Meters
	double v_y = 0.25;											// m/s
	double depletionTime = 200.;							// Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 9, yCells = 27;
	double xc, yc, s;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string solPath = getDataPath() + "caseStudy/msrLumpDepletion/";
	std::string outputFileName = "msr2DDepletionMedium9x27.out";
	std::string outputFileNameMatrix = solverType+"MSR2DDepletionMedium9x27.csv";
	std::string limiter = "First order upwind";
	std::ofstream outputFile;
	outputFile.open(outputFileName, std::ios_base::app);
	MatrixD refSolData;
	FILE * pOutputFile;
	meshCell* cell = nullptr;
	std::vector<std::string> wallTransNames={"Nb-93", "Mo-95", "Mo-97", "Mo-98", "Mo-99",
													     "Mo-100", "Tc-99", "Ru-101", "Ru-102", "Ru-103", 
														  "Ru-104", "Ru-106", "Rh-105", "Rh-103", "Pd-105",
														  "Pd-107", "Pd-108", "Ag-109", "Sn-126"};
	std::vector<std::string> gasNames = {"H-1Gas", "N-14Gas", "O-16Gas", "Kr-83Gas", 
													 "Xe-131Gas", "Xe-133Gas", "Xe-135Gas"};
	// {graphite, heat exchanger, piping}
	std::vector<double> massTransCoeffs = {5.334e-6, 4.657e-5, 1.041e-4};

	if (myid == 0){
		pOutputFile = fopen(outputFileName.c_str(), "a");
		fprintf(pOutputFile, "Solver: %s \n", solverType.c_str());
		fprintf(pOutputFile, "Total problem time: %4.3e\n", totalTime);
		fprintf(pOutputFile, "yLength: %4.3e\n", yLength);
		fprintf(pOutputFile, "xLength: %4.3e\n", xLength);
		fprintf(pOutputFile, "dt: %8.7e\n", dt);
		fprintf(pOutputFile, "yCells: %d\n", yCells);
		fprintf(pOutputFile, "xCells: %d\n", xCells);
		fclose(pOutputFile);
	}

	// Sets file paths for the input data	
	std::string speciesNamesFile = path + "speciesInputNamesMSRMassTransportMedium.dat";
	std::string speciesDecayFile = path + "speciesInputDecayMSRMassTransportMedium.dat";
	std::string speciesTransFile = path + "speciesInputTransMSRMassTransportMedium.dat";

	// Builds the model
	modelMesh model(xCells, yCells, xLength, yLength);

	// Adds boundary surface
	model.addBoundarySurface("east");
	model.addBoundarySurface("west");

	// set y velocity
	model.setConstantYVelocity(v_y);

	// Builds the species object
	speciesDriver spec = speciesDriver(&model);

	// Sets the variables
	model.setSystemTemperature(path + "temperatureMap9x27.csv");
	model.setSystemNeutronFlux(path + "neutronFluxMap9x27.csv");
	model.setSystemGasVoidFraction(path + "voidMap9x27.csv");
	model.setSystemGasInterfacialAreaCon(path + "gasIntAreaMap9x27.csv");
	model.setSystemWallInterfacialAreaCon(path + "wallIntAreaMap9x27.csv");

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);

	// Sets the flux limiter type
	spec.setFluxLimiter(limiter);

	// Adds the speices
	ids = spec.addSpeciesFromFile(speciesNamesFile);
	spec.setGasSpargingFromFile(path + "speciesInputGasSpargingMedium.txt");

	// Loops through to build the liqIDs and wallIDs for wall deposition 
	std::vector<int> liqIDs = {}, wallIDs = {};
	std::vector<bool> infSinks = {};
	for (int k = 0; k < wallTransNames.size(); k++){
		std::string name = wallTransNames[k];
		liqIDs.push_back(spec.getSpeciesID(name+"Liq"));
		wallIDs.push_back(spec.getSpeciesID(name+"Wall"));
		infSinks.push_back(true);
	}

	// Loops through the model to set the wall mass transport coefficients
	for (int i = 0; i < xCells; i++){
		for (int j = 0; j < yCells; j++){
			// The core region. This contains graphite
			if (j > 0 and j < 9){
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[0]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}
			// Heat exchanger
			else if (j >= 12 and j <= 14){
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[1]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}
			else{
				std::vector coeffs(wallTransNames.size(), massTransCoeffs[2]);
				spec.setWallDeposition(i, j, coeffs, liqIDs, wallIDs, infSinks);
			}

			// Adds the gas removal
			if (j == 10){
				for (int k = 0; k < gasNames.size(); k++){
					spec.setRemoval(i, j, spec.getSpeciesID(gasNames[k]), -2.05e-4);
				}
			}
		}
	}
	
	// Sets the boundary condition values
	std::vector<double> bcs = {};

	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	//readCSV(refSolData, std::string(path + "solutionMSR2DDepletionSmall.csv"));

	// Print species name to file
	if (myid == 0){
		pOutputFile = fopen(outputFileName.c_str(), "a");
		fprintf(pOutputFile, "%s %s %s ", "variables", "x", "y");
		for (int id; id < ids.size(); id++){
			std::string name = spec.getSpeciesName(id);
			fprintf(pOutputFile, "%s ", name.c_str());
		}
		fprintf(pOutputFile, "\n");
		fclose(pOutputFile);
	}

	// Sets the species sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);

	// Adds the boundary conditions
	spec.setBoundaryCondition("newmann","east", ids, bcs);
	spec.setBoundaryCondition("newmann","west", ids, bcs);
	spec.setBoundaryCondition("periodic","south", ids);
	spec.setBoundaryCondition("periodic","north", ids);

	// Writes transition matrix and initial condition
	if (myid == 0){
		spec.writeTransitionMatrixToFile("transitionMatrixMSR2DDepletionMedium9x27.csv");
		spec.writeInitialConditionToFile("initialConditionMSR2DDepletionMedium9x27.csv");
		pOutputFile = fopen(outputFileName.c_str(), "a");
	}

	// Loops over the time steps to solve
	double totalSolveTime = 0;
	for (int k = 0; k < steps; k++){
		auto start = std::chrono::high_resolution_clock::now();
		t = t + dt;
		//std::cout << "solving" << std::endl;
		spec.solve(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		if (myid == 0){
			solveTime = duration.count()/1.e6;
			totalSolveTime += solveTime;

			fprintf(pOutputFile, "time: %5.4e\n", t);
			int index = 0;
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					fprintf(pOutputFile, "%5.4e %5.4e ", xc, yc);
					cell = model.getCellByLoc(i,j); xc = cell->x; yc = cell->y;

					for (int id = 0; id < ids.size(); id++){
						std::string name = spec.getSpeciesName(ids[id]);
						double con = spec.getSpecies(i, j, ids[id]);
						fprintf(pOutputFile, "%17.16e ", con);
						printf("%5.4e %5.4e %s %17.16e \n", xc, yc, name.c_str(), con);
						solData(index, k) = con;
						index ++;
					}
					fprintf(pOutputFile, "\n");
				}
			}
			//VectorD refSolVect = refSolData.col(k), solVect = solData.col(k);
			//double E1 = computeRelativeE1(refSol, sol);
			//double E2 = computeRelativeE2(refSol, sol);
			//double Einf = computeRelativeEinfty(refSol, sol);
			double E1 = 0.;
			double E2 = 0.;
			double Einf  = 0.;
			printf("%s %e %e %e %e\n", solverType.c_str(), t, Einf, E1, E2);
		}
	}
	if (myid == 0){
		printf("%s Solve Time: %f\n", solverType.c_str(), totalSolveTime);
		fprintf(pOutputFile, "end\n");
		writeCSV(solData, outputFileNameMatrix);
	}
}
//*****************************************************************************
// Main test
//*****************************************************************************
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;
	std::vector<std::string> solvers {"CRAM", "hyperbolic", "parabolic",
	"pade-method1", "pade-method2", "taylor"};
	//std::vector<std::string> solvers {"BDF1", "BDF2", "BDF3", "BDF4", "BDF5", "BDF6"};
	//std::vector<std::string> solvers {"BDF1"};
	//std::vector<std::string> solvers {"forward euler", "explicit midpoint", "kutta third-order", 
	//"classic fourth-order"};
	//std::vector<std::string> solvers {"explicit midpoint", "forward euler"};
	//std::vector<std::string> solvers {"CRAM"};
	//std::vector<std::string> solvers {"CRAM", "hyperbolic", "parabolic"};
	//std::vector<std::string> solvers {"hyperbolic", "parabolic"};
	//std::vector<std::string> solvers {"pade-method1"};

	// Loops over different solvers
	for (std::string &solverType : solvers){
		//if(myid == 0){std::cout << solverType << std::endl;};
		//msrLumpDepletion(myid, solverType);

		// For my M&C paper this does not include mass transport
		//msr2DDepletion(myid, solverType);

		// Neutron precursor problem from my dissertation
		//neutronPrecursors(myid, solverType);
		// These are for my dissertation these include mass transport
		//msrDepletionSmallLumped(myid, solverType);
		//msrDepletionMediumLumped(myid, solverType);
		msr2DDepletionSmall3x9(myid, solverType);
		//msr2DDepletionMedium3x9(myid, solverType);
		//msr2DDepletionSmall9x27(myid, solverType);
		//msr2DDepletionMedium9x27(myid, solverType);
	}

	mpi.finalize();
}
