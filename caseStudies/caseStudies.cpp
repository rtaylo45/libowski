#include "mpiProcess.h"
#include "modelMesh.h"
#include "speciesDriver.h"
#include "meshCell.h"
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
	double depletionTime = 500.; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps;
	int xCells = 1, yCells = 1;
	double xLength = 220.0/1000., yLength = 220.0/1000.;
	std::vector<int> ids;
	std::string path = getDataPath() + "caseStudy/";
	std::string outputFileName = "caseStudy1.out";
	//std::string outputFileNameMatlab = "caseStudy1"+solverType+"Substeps12.csv";
	std::string outputFileNameMatlab = "caseStudy1"+solverType+".csv";
	std::ofstream outputFile, outputFileMatlab;
	outputFile.open(outputFileName, std::ios_base::app);
	MatrixD refSolData;


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

	// Get solution
	readCSV(refSolData, std::string(path+"caseStudy1Solution.csv"));

	outputFile << "solverName: " << solverType << std::endl;
	outputFile.precision(16); 

	// Writes transition matrix and initial condition
	spec.writeTransitionMatrixToFile("transitionMatrixCaseStudy1.csv");
	spec.writeInitialConditionToFile("initialConditionCaseStudy1.csv");
	MatrixD solData = MatrixD::Zero(ids.size()+1,10);

	// Loops over the time steps to solve
	for (int k = 0; k < steps; k++){
		auto start = std::chrono::high_resolution_clock::now();
		t = t + dt;
		spec.solve(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		if (myid == 0){
			outputFile << "Time: " << t << std::endl;
			std::cout << "Time: " << t << " Solve Time: "
			<< duration.count()/1.e6 << std::endl;
			for (int id = 0; id < ids.size(); id++){
				std::string name = spec.getSpeciesName(ids[id]);
				double con = spec.getSpecies(0, 0, ids[id]);
				outputFile << name << " " << con << std::endl;
				solData(id, k) = con;
				double error = std::abs(refSolData(id, k) - con)/refSolData(id,k);
				//std::cout << k << " " << id << " " << error << std::endl;
			}
		}
	}
	if (myid == 0){
		double rmse = computeRelativeRMSE(refSolData, solData);
		assert(rmse < 1.e-7);
		printf("%s RMSE %e\n", solverType.c_str(), rmse);
	}
	outputFile << " " << std::endl;
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
	std::string outputFileName = "caseStudy2.out";
	std::string path = getDataPath() + "caseStudy/";
	std::ofstream outputFile;
	VectorD refSolData;
	outputFile.open(outputFileName, std::ios_base::app);

	// Files for the source terms and species names	
	std::string speciesNamesFile = path + "speciesInputNames.dat";
	std::string speciesDecayFile = path + "speciesInputDecay.dat";
	std::string speciesTransFile = path + "speciesInputTrans.dat";

	// Builds the model mesh
	modelMesh model(xCells, yCells, xLength, yLength);

	// Inits the species driver
	speciesDriver spec = speciesDriver(&model);

	// Species IDs
	ids = spec.addSpeciesFromFile(speciesNamesFile);

	// Sets all the decay and transmutation sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);
	VectorD solData = VectorD::Zero(xCells*yCells*ids.size()+1);
	readCSV(refSolData, std::string(path+"caseStudy2Solution.csv"));

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

	// Adds periodic boundary condions to all isotopes
	spec.setBoundaryCondition("periodic","south", ids);
	spec.setBoundaryCondition("periodic","north", ids);

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);

	if (myid == 0){
		outputFile << "solverName: " << solverType << std::endl;
		outputFile.precision(16); 
		spec.writeTransitionMatrixToFile("transitionMatrixCaseStudy2.csv");
		spec.writeInitialConditionToFile("initialConditionCaseStudy2.csv");
	}


	// Loops to solve the problem
	for (int k = 0; k < steps; k++){
		t = t + dt;
		if (myid == 0){std::cout << t << std::endl;};
		spec.solve(t);
	}
	if (myid == 0){
		outputFile << "Time: " << t << std::endl;
		outputFile << "Variables: " << "i " << "j " << "SpecId " << 
			"SpecName " << "Con" << std::endl;
		// Loops to print results
		int index = 0;
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				// Loops over the species 
				for (int id = 0; id < ids.size(); id++){
					std::string name = spec.getSpeciesName(ids[id]);
					double con = spec.getSpecies(i, j, ids[id]);
					outputFile << i << " " << j << " " << id << " " << name 
						<< " " << con << std::endl;
					solData[index] = con;
					index ++;
				}
			}
		}
		outputFile << " " << std::endl;
		if (myid == 0){
			double rmse = computeRelativeRMSE(refSolData, solData);
			printf("%s RMSE %e\n", solverType.c_str(), rmse);
		}
	}
}

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
	int steps = 6;
	double totalTime = 60.0;
	double dt = totalTime/steps;
	//std::vector<double> timeSteps = {0.1, 0.5, 1.0, 10.0, 20.0};
	//int xCells = 5, yCells = 20;
	int xCells = 5, yCells = 20;
	double xLength = 50.0, yLength = 400.0;
	double v_y = 25.0, v_x = 0.0;
	//double v_y = 60.0, v_x = 0.0;
	double xc, yc, s, g;
	MatrixD refSolData;
	std::vector<int> ids;
	std::vector<double> decay = {0.0127, 0.0317, 0.115, 0.311, 1.4, 3.87};
	std::vector<double> beta = {0.06, 0.364, 0.349, 0.628, 0.179, 0.07};
	std::vector<double> bcs = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	std::string path = getDataPath() + "caseStudy/";
	std::string speciesNamesFile = path + "neutronPrecursorsInputNames.dat";
	std::string limiter = "First order upwind";
	std::string outputFileName = "caseStudyNeutronPrecursors.out";
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

	// Sets the flux limiter type
	spec.setFluxLimiter(limiter);

	// Species IDs
	ids = spec.addSpeciesFromFile(speciesNamesFile);

	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	readCSV(refSolData, std::string(path +
		"caseStudyNeutronprecursorsSolution.csv"));

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
		//t = timeSteps[k];
		auto start = std::chrono::high_resolution_clock::now();
		spec.solve(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		//std::cout << "Time: " << t << " Solve Time: "
		//<< duration.count()/1.e6 << std::endl;
		totalSolvetime += duration.count()/1.e6;

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
			printf("%s dt = %f RMSE %e\n", solverType.c_str(), t, rmse);
			if (solverType != "pade-method1"){
				assert(rmse < 1.e-10);
			}
			else{
				assert(rmse < 1.e-2);
			}
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
	int steps = 10;
	double depletionTime = 500.; // Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 1, yCells = 1;
	double xLength = 1., yLength = 1., xc, yc;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string outputFileName = "msrMediumLumpDepletion.out";
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
	std::string speciesNamesFile = path + "speciesInputNamesMSRMassTransportMedium.dat";
	std::string speciesDecayFile = path + "speciesInputDecayMSRMassTransportMedium.dat";
	std::string speciesTransFile = path + "speciesInputTransMSRMassTransportMedium.dat";

	// Builds the model
	modelMesh model(xCells, yCells, xLength, yLength);

	// Builds the species object
	speciesDriver spec = speciesDriver(&model);

	// Sets the neutron flux
	model.setSystemNeutronFlux(1.e13);
	model.setSystemGasInterfacialAreaCon(150.0);
	model.setSystemTemperature(900.0);
	model.setSystemGasVoidFraction(0.003);
	model.setSystemWallInterfacialAreaCon(100.0);

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);

	// Adds the speices
	ids = spec.addSpeciesFromFile(speciesNamesFile);

	// Sets the species sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);
	spec.setGasSpargingFromFile(path + "speciesInputGasSpargingMedium.txt");
	spec.setWallDepositionFromFile(path + "speciesInputWallDepositionMedium.txt");

	// Writes transition matrix and initial condition
	//spec.writeTransitionMatrixToFile("transitionMatrixMSRMediumLumpDepletion.csv");
	//spec.writeInitialConditionToFile("initialConditionMSRMediumLumpDepletion.csv");
	pOutputFile = fopen(outputFileName.c_str(), "a");

	// Gets the solution
	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	readCSV(refSolData, std::string(path + "msrMediumLumpDepletionSolution.csv"));

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

			cell = model.getCellByLoc(0,0); xc = cell->x; yc = cell->y;

			fprintf(pOutputFile, "time: %5.4e\n", t);

			for (int id = 0; id < ids.size(); id++){
				std::string name = spec.getSpeciesName(ids[id]);
				double con = spec.getSpecies(0, 0, ids[id]);
				fprintf(pOutputFile, "%5.4e %5.4e %s %17.16e \n", xc, yc, name.c_str(), con);
				printf("%5.4e %5.4e %s %17.16e \n", xc, yc, name.c_str(), con);
				solData(id, k) = con;
			}
			fprintf(pOutputFile, "\n");
			//VectorD refSolVect = refSolData.col(k), solVect = solData.col(k);
			//double rmse = computeRelativeRMSE(refSolVect, solVect);
			double rmse = 0.0;
			printf("Solve Step: %d %s Solve Time: %f RMSE %e\n", k, solverType.c_str(), solveTime, rmse);
		}
	}
	if (myid == 0){
		fprintf(pOutputFile, "end\n");
		//double rmse = computeRelativeRMSE(refSolData, solData);
		double rmse = 0.0;
		printf("%s Solve Time: %f RMSE %e\n", solverType.c_str(), totalSolveTime, rmse);
	}
}

//*****************************************************************************
// MSR pipe depletion
//*****************************************************************************
void msrPipeDepletion(int myid, std::string solverType){
	double t = 0.0;
	int steps = 5;
	double xLength = 0.0508, yLength = 1.7018;		// Meters
	double v_y = 0.25;											// m/s
	double depletionTime = 50.;							// Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 1, yCells = 10;
	double xc, yc;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string outputFileName = "msrMediumPipeDepletion.out";
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
	std::string speciesNamesFile = path + "speciesInputNamesMSRMassTransportMedium.dat";
	std::string speciesDecayFile = path + "speciesInputDecayMSRMassTransportMedium.dat";
	std::string speciesTransFile = path + "speciesInputTransMSRMassTransportMedium.dat";

	// Builds the model
	modelMesh model(xCells, yCells, xLength, yLength);

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
	model.setConstantYVelocity(v_y);


	// Builds the species object
	speciesDriver spec = speciesDriver(&model);

	// Sets the matrix exp solver
	spec.setMatrixExpSolver(solverType);

	// Adds the speices
	ids = spec.addSpeciesFromFile(speciesNamesFile);

	// Sets the species sources
	spec.setSpeciesSourceFromFile(speciesDecayFile, speciesTransFile);

	// Adds periodic boundary condions to all isotopes
	spec.setBoundaryCondition("periodic","south", ids);
	spec.setBoundaryCondition("periodic","north", ids);

	// Writes transition matrix and initial condition
	spec.writeTransitionMatrixToFile("transitionMatrixMSRMediumPipeDepletion.csv");
	spec.writeInitialConditionToFile("initialConditionMSRMediumPipeDepletion.csv");
	pOutputFile = fopen(outputFileName.c_str(), "a");

	// Gets the solution
	MatrixD solData = MatrixD::Zero(xCells*yCells*ids.size()+1, steps);
	//readCSV(refSolData, std::string(path + "msrPipeDepletionSolution.csv"));

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

			cell = model.getCellByLoc(0,0); xc = cell->x; yc = cell->y;

			fprintf(pOutputFile, "time: %5.4e\n", t);

			for (int id = 0; id < ids.size(); id++){
				std::string name = spec.getSpeciesName(ids[id]);
				double con = spec.getSpecies(0, 0, ids[id]);
				fprintf(pOutputFile, "%5.4e %5.4e %s %17.16e \n", xc, yc, name.c_str(), con);
				solData(id, k) = con;
			}
			fprintf(pOutputFile, "\n");
			//VectorD refSolVect = refSolData.col(k), solVect = solData.col(k);
			//double rmse = computeRelativeRMSE(refSolVect, solVect);
			double rmse = 0.0;
			printf("Solve Step: %d %s Solve Time: %f RMSE %e\n", k, solverType.c_str(), solveTime, rmse);
		}
	}
	if (myid == 0){
		fprintf(pOutputFile, "end\n");
		//double rmse = computeRelativeRMSE(refSolData, solData);
		double rmse = 0.0;
		printf("%s Solve Time: %f RMSE %e\n", solverType.c_str(), totalSolveTime, rmse);
	}
}

//*****************************************************************************
// MSR 2-D depletion
//*****************************************************************************
void msr2DDepletion(int myid, std::string solverType){
	double t = 0.0;
	int steps = 1;
	double coreLength = 1.7018;
	double xLength = 0.0508, yLength = 4.*coreLength;		// Meters
	double v_y = 0.25;											// m/s
	double depletionTime = 1.;							// Days
	double totalTime = depletionTime*24.*60.*60.;
	double dt = totalTime/steps, solveTime;
	int xCells = 10, yCells = 20;
	double xc, yc, s;
	std::vector<int> ids;
	std::string path = getDataPath() + "msr/";
	std::string outputFileName = "msr2DDepletionSmall.out";
	std::string limiter = "First order upwind";
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
				s = 1.e13*sy*sx;
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
		spec.writeTransitionMatrixToFile("transitionMatrixMSRPipeDepletion.csv");
		spec.writeInitialConditionToFile("initialConditionMSRPipeDepletion.csv");
		pOutputFile = fopen(outputFileName.c_str(), "a");
	}

	// Loops over the time steps to solve
	double totalSolveTime = 0;
	for (int k = 0; k < steps; k++){
		auto start = std::chrono::high_resolution_clock::now();
		t = t + dt;
		std::cout << "solving" << std::endl;
		spec.solve(t);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end - start);
		if (myid == 0){
			solveTime = duration.count()/1.e6;
			totalSolveTime += solveTime;

			cell = model.getCellByLoc(0,0); xc = cell->x; yc = cell->y;

			fprintf(pOutputFile, "time: %5.4e\n", t);
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					fprintf(pOutputFile, "%5.4e %5.4e ", xc, yc);

					for (int id = 0; id < ids.size(); id++){
						std::string name = spec.getSpeciesName(ids[id]);
						double con = spec.getSpecies(0, 0, ids[id]);
						fprintf(pOutputFile, "%17.16e ", con);
					}
					fprintf(pOutputFile, "\n");
				}
			}
		}
	}
	if (myid == 0){
		printf("%s Solve Time: %f\n", solverType.c_str(), totalSolveTime);
		fprintf(pOutputFile, "end\n");
	}
}

//*****************************************************************************
// Main test
//*****************************************************************************
int main(){
	int myid = mpi.rank;
	int numprocs = mpi.size;
	//std::vector<std::string> solvers {"CRAM", "hyperbolic", "parabolic",
	//"pade-method2", "taylor"};
	std::vector<std::string> solvers {"hyperbolic"};
	//std::vector<std::string> solvers {"pade-method1", "pade-method2"};

	// Loops over different solvers
	for (std::string &solverType : solvers){
		//if(myid == 0){std::cout << solverType << std::endl;};
		singleCellDepletion(myid, solverType);
		//if (solverType != "taylor"){
		//	pipeDepletion(myid, solverType);
			//msrPipeDepletion(myid, solverType);
		//}
		//neutronPrecursors(myid, solverType);
		msrLumpDepletion(myid, solverType);
		//msr2DDepletion(myid, solverType);
	}

	mpi.finalize();
}
