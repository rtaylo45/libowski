#include <iostream>
#include "sys.h"
#include "exception.h"
#include "parser.h"
#include "modelMesh.h"
#include "speciesDriver.h"
using namespace std;

int main(int argc, char *argv[]){

	// Error checking for passed arguments
	if (argc == 0){
		string errorMessage =
			" No input file\n";
		libowskiException::runtimeError(errorMessage);
	}
	else if (argc > 2){
		string errorMessage =
			" Too many arguments passed\n";
		libowskiException::runtimeError(errorMessage);
	}

	// Gets the file name
	string inputFile = argv[1];
	parser fileParser;
  // Parses the input file and creates the data blocks
	fileParser.parseFile(inputFile);
  // Parses the mesh block
  modelMesh* meshPtr = fileParser.parseMeshBlock();
  // Parses the species block
  speciesDriver* speciesdriverPtr = fileParser.parseSpeciesBlock(meshPtr);
  // Parses the auxvariable block
  fileParser.parseAuxVariableBlock(meshPtr);
	// Creats the problem object
	// Problem problem;
	// Parses the file and builds the objects
	// problem.parseFile(fileName)
	// Runs the problem
	// problem.run()
	// Adds output for the problem
	// problem.output()
	// Cleans up the problem
	// problem.clean()

	return 0;
}
