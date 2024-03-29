//*****************************************************************************
// Author: Zack Taylor
//
// The model mesh. Houses info on the problem domain
//*****************************************************************************
#ifndef MODELMESH_H
#define MODELMESH_H
#include "meshCell.h"
#include "meshCellFace.h"
#include "matrixTypes.h"
#include "utilBase.h"
#include "surface.h"
#include <vector>
#include <assert.h>
#include <string>

class modelMesh {

  public:
  // Total length in x direction [m]
  double xLength = 0.0;
  // Total length in y direction [m]
  double yLength = 0.0;
  // Number of cells in the x direction
  int numOfxCells = 0;
  // Number of cells in the y direction
  int numOfyCells = 0;
  // Total number of cells
  int numOfTotalCells = 0;
  // Vector of all cells
  std::vector<meshCell> meshCells;
  // Vector of all cell faces
  std::vector<meshCellFace> meshCellFaces;
  // Vector of cell surfaces
  std::vector<surface> surfaces;
  // Change in x direction [m]
  double dx = 0.0;
  // Change in y direction [m]
  double dy = 0.0;

  public:
  // Constructor
  modelMesh(int, int, double, double);
  // Builds the geometry
  void buildGeometry();
  // Gets the node by location from i,j
  meshCell* getCellByLoc(int, int);
  // Gets the node by location from absolution index k
  meshCell* getCellByLoc(int);
  // Sets a constant x velocity across the whole problem
  void setConstantXVelocity(double);
  // Sets a constant x velocity across a column of cells
  void setConstantXVelocity(double, int);
  // Sets a constant y velocity across the whole problem
  void setConstantYVelocity(double);
  // Sets a constant y velocity across a row of cells
  void setConstantYVelocity(double, int);
  // Sets a temperature in the whole system
  void setSystemTemperature(double);
  // Sets a temperature in the whole system from a file
  void setSystemTemperature(const std::string);
  // Sets a pressure in the whole system
  void setSystemPressure(double);
  // Sets a pressure in the whole system from a file
  void setSystemPressure(const std::string);
  // Sets a neutron flux in the whole system
  void setSystemNeutronFlux(double);
  // Sets a neutron flux in the whole system from a file
  void setSystemNeutronFlux(const std::string);
  // Sets the gas interfacial area concentation in the whole system
  void setSystemGasInterfacialAreaCon(double);
  // Sets the gas interfacial area concentation in the whole system from a file
  void setSystemGasInterfacialAreaCon(const std::string);
  // Sets the wall interfacial area concentation in the whole system
  void setSystemWallInterfacialAreaCon(double);
  // Sets the wall interfacial area concentation in the whole system from a file
  void setSystemWallInterfacialAreaCon(const std::string);
  // Sets the gas void fraction in the whole system
  void setSystemGasVoidFraction(double);
  // Sets the gas void fraction in the whole system from a file
  void setSystemGasVoidFraction(const std::string);
  // Set temperature in a cell
  void setCellTemperature(int, int, double);
  // Set pressure in a cell
  void setCellPressure(int, int, double);
  // Set neutron flux in a cell
  void setCellNeutronFlux(int, int, double);
  // Sets the gas interfacial area concentation in a cell
  void setCellGasInterfacialAreaCon(int, int, double);
  // Sets the wall interfacial area concentation in a cell
  void setCellWallInterfacialAreaCon(int, int, double);
  // Sets the gas void fraction in a cell
  void setGasVoidFraction(int, int, double);
  // Adds a physical surface to a cell
  void addSurface(int, int, std::string);
  // Adds a physical surface along a boundary
  void addBoundarySurface(std::string);
  // Cleans the model
  void clean();

  private:
  // Creates the cells
  void createCells();
  // Connects the cells
  void connectCells();
  // Creates the cell faces
  void createCellFaces();
  // Checks the i,j for validity
  bool checkCellLoc(int, int);

};
#endif
