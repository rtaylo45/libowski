//*****************************************************************************
// Author: Zack Taylor
// 
//*****************************************************************************
#ifndef PARSER_H
#define PARSER_H

#include <memory>
#include "modelMesh.h"
#include "speciesDriver.h"

struct dataBlock {

  std::map<std::string, std::vector<std::string>> variableValueMap;
  std::string blockName;

  dataBlock(std::string _blockName){
    blockName = _blockName;
  }
  void addVariable(std::string varName, std::vector<std::string> varVal){
    variableValueMap[varName] = varVal;
  }
  std::vector<std::string> getVariableNames(){
    std::vector<std::string> vars;
    for (auto const& element : variableValueMap){
      vars.push_back(element.first);
    }
    return vars;
  }
  std::vector<std::string> getVariableValues(std::string varName){
    return variableValueMap[varName];
  }
};

class parser {

  private:
  std::map<std::string, dataBlock> inputDeckBlocks;

  public:
  // Constructor
  parser(){};
  // Reads the file and builds the input blocks
  void parseFile(std::string); 
  // Gets the mesh model for the problem
  std::unique_ptr<modelMesh> getModelMesh();
  // Gets the species driver for the problem
  //unique_ptr<speciesDriver> getSpeciesDriver();

  private:
  // Parses the mesh block and creates a pointer
  std::unique_ptr<modelMesh> parseMeshBlock();
  // Parse Mesh variable block and adds it to the model mesh pointer
  //void parseMeshVariableBlock(unique_ptr<modelMesh>&);
  // Parse species block and creates a pointer
  //unique_ptr<speciesDriver> parseSpeciesBlock();
};

#endif
