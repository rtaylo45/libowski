#include <iostream>
#include "sys.h"
#include "exception.h"
#include "parser.h"

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
	FileParser parser;
	parser.readFile(inputFile);
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
