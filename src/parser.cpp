//*****************************************************************************
// Author: Zack Taylor
//*****************************************************************************
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include "parser.h"
#include "sys.h"
#include "utilBase.h"
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
  string blockName;
	while (inFile){

		string line;
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
          if (blockName == "Mesh"){
            string delimiter = "=";
            size_t pos = 0;
            // Get rid of white space
            line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
            //vector<string> splitLine splitStr(line, delimiter);
            //for (auto i : splitLine) cout << i << endl;
          }  
				}
			}
		}
	}
	
	inFile.close();	
}














