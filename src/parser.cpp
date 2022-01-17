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
#include "exception.h"

using namespace std;

//**************************************************************************
//
// @param blockName Name of the data block
//**************************************************************************
dataBlock* parser::getDataBlock(const std::string &blockName ){
  try {
    return &inputDeckBlocks.at(blockName);
  }
  catch (const std::out_of_range& e) {
    dataBlock newDat = dataBlock(blockName);
    inputDeckBlocks.insert({blockName, newDat});
    return &inputDeckBlocks.at(blockName);

  }
}

//**************************************************************************
// Parses a file
//
// @param fname Input file to read
//**************************************************************************
void parser::parseFile(const std::string &fname){

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
          if (anyIn(blockName, blockNames)){
            dataBlock* datPtr = getDataBlock(blockName);
            // Get rid of white space
            line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
            vector<string> splitLine = splitStr(line, string("[:,=]"));
            string varName = string(splitLine[0]);
            splitLine.erase(splitLine.begin());
            datPtr->addVariable(varName, splitLine);
          }
          else {
            string errorMessage =
              " The input block you have passed is not recognized by libowski.";
              libowskiException::runtimeError(errorMessage);
          }
        }
      }
    }
  }
  inFile.close();
}

//**************************************************************************
// Parses the variables in the mesh block and returns a pointer to the
// created mesh variable.
//**************************************************************************
modelMesh* parser::parseMeshBlock(){
  dataBlock* datPtr = getDataBlock(string("Mesh"));
  // Gets the data
  double xLength = stod(datPtr->getVariableValues(string("xLength"))[0]);
  double yLength = stod(datPtr->getVariableValues(string("yLength"))[0]);
  int xCells = stoi(datPtr->getVariableValues(string("xCells"))[0]);
  int yCells = stoi(datPtr->getVariableValues(string("yCells"))[0]);
  // Creates the pointer
  modelMesh* meshPtr = new modelMesh(xCells, yCells, xLength, yLength);
  return meshPtr;
}

//**************************************************************************
// Parses the variables in the species block and returns a pointer to the
// created speciesDriver variable.
//**************************************************************************
speciesDriver* parser::parseSpeciesBlock(modelMesh* modelMeshPtr){
  dataBlock* datPtr = getDataBlock(string("Species"));
  // Gets data
  vector<string> speciesNames = datPtr->getVariableValues(string("name"));
  vector<string> speciesPhases = datPtr->getVariableValues(string("phase"));
  vector<string> speciesDiffusivity = datPtr->getVariableValues(string("diffusivity"));
  vector<string> speciesMolarMasses = datPtr->getVariableValues(string("molar_mass"));
  vector<string> speciesInitialCon = datPtr->getVariableValues(string("initial_condition"));

  // Does some checks on the variables
  if (speciesNames.size() != speciesPhases.size())
    libowskiException::runtimeError("Species names and phases are not the same size");
  if (speciesNames.size() != speciesMolarMasses.size())
    libowskiException::runtimeError("Species names and molar masses are not the same size");
  if (speciesDiffusivity.size() != 0 ){
    if (speciesDiffusivity.size() != speciesNames.size())
      libowskiException::runtimeError("Species names and diffusivity are not the same size");
  }
  if (speciesInitialCon.size() != 0 ){
    if (speciesInitialCon.size() != speciesNames.size())
      libowskiException::runtimeError("Species names and initial condition are not the same size");
  }

  // Creates the pointer
  speciesDriver* speciesDriverPtr = new speciesDriver(modelMeshPtr);

  // Adds the species to the problem
  for (int i = 0; i < speciesNames.size(); i++){
    string spec_name = speciesNames[i] + "_" + speciesPhases[i];
    double molarMass = stod(speciesMolarMasses[i]);
    double diffusivity = 0.0;
    if (speciesDiffusivity.size())
      diffusivity = stod(speciesDiffusivity[i]);
    double initialCondition = 0.0;
    if (speciesInitialCon.size())
      initialCondition = stod(speciesInitialCon[i]);
    bool transported = true;
    if (speciesPhases[i] == "solid")
      transported = false;
    int specID = speciesDriverPtr->addSpecies(molarMass, ititialCondition,
      diffusivity, spec_name, transported);
  }
  return speciesDriverPtr;
}
