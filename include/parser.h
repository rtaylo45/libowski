//*****************************************************************************
// Author: Zack Taylor
//
//*****************************************************************************
#ifndef PARSER_H
#define PARSER_H

#include <memory>
#include <string>
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
  void print(){
    for (auto const & element : variableValueMap){
      printf("Variable Name: %s Varaible Values: ", element.first.c_str());
      for (auto const & var : element.second){
        printf("%s ", var.c_str());
      }
      printf("\n");
    }
    printf("\n");
  }
};

class parser {

  private:
  std::map<std::string, dataBlock> inputDeckBlocks;

  public:
  // Constructor
  parser(){};
  // Reads the file and builds the input blocks
  void parseFile(const std::string&);
  // Gets the data block from its name
  dataBlock* getDataBlock(const std::string&);
  // Gets the mesh model for the problem
  std::unique_ptr<modelMesh> getModelMesh();
  // Gets the species driver for the problem
  //unique_ptr<speciesDriver> getSpeciesDriver();

  private:
  // List of the accepted block names
  const std::vector<std::string> blockNames = {"Mesh", "Species"};
  // Parses the mesh block and creates a pointer
  std::unique_ptr<modelMesh> parseMeshBlock();
  // Parse Mesh variable block and adds it to the model mesh pointer
  //void parseMeshVariableBlock(unique_ptr<modelMesh>&);
  // Parse species block and creates a pointer
  //unique_ptr<speciesDriver> parseSpeciesBlock();
};

#endif
