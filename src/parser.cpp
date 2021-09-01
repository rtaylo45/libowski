//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include "parser.h"
#include "sys.h"
using namespace std;

//**************************************************************************
// Parses a file
//
// @param fname Input file to read
//**************************************************************************
void parser::parseFile(std::string fname){

	// Checks to see if the file exist
	checkFileExists(fname);

	ifstream inFile(fname);
	bool inBlock = false;
	while (inFile){

		string line, blockName;
		getline(inFile, line);
		if (line.size() > 0){
			// This is a comment line
			if (line.at(0) == '#'){
				continue;
			}
			// Start of a block
			else if (line.at(0) == '['){
				line.erase(remove(line.begin(), line.end(), '['), line.end());	
				line.erase(remove(line.begin(), line.end(), ']'), line.end());
				blockName = line;
				inBlock = true;
			}
			// end of a block
			else if (line.at(1) == ']'){
				inBlock = false;
			}
			else {
				if (inBlock){
					cout << blockName << " " << line << endl;

				}
			}
		}
	}
	
	inFile.close();	
}














