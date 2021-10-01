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
//
// @param blockName Name of the data block
//**************************************************************************
dataBlock parser::getDataBlock(const & std::string blockName){

}

//**************************************************************************
// Parses a file
//
// @param fname Input file to read
//**************************************************************************
void parser::parseFile(const & std::string fname){

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
          // Parses info specific to the mesh block
          if (blockName == "Mesh"){
            dataBlock dat = getDataBlock(blockName);
            string delimiter = "[=]";
            size_t pos = 0;
            // Get rid of white space
            line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
            vector<string> splitLine = splitStr(line, delimiter);
            for (auto i : splitLine) {
              dat.addVariable(splitLine[0], {splitLine[1]});
            }
            dat.print();
            //inputDeckBlocks[blockName] = dat;
          }
          else if (blockName == "Species"){
            dataBlock dat = getDataBlock(blockName);
            size_t pos = 0;
            // Get rid of white space
            line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
            vector<string> splitLine = splitStr(line, string("[=,]"));
            string varName = string(splitLine[0]);
            splitLine.erase(splitLine.begin());
            dat.addVariable(varName, splitLine);
            dat.print();
            //inputDeckBlocks[blockName] = dat;
          }  
				}
			}
		}
	}
	inFile.close();	
}














