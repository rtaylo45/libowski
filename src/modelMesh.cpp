//*****************************************************************************
// Author: Zack Taylor

// The mesh class. Houses info on the problem domain
//*****************************************************************************
#include "modelMesh.h"

//**************************************************************************
// Constructor
//
// @param xCells    Number of cells in the x direction
// @param yCells    Number of cells in the y direction
// @param xDirLength  Total length in the x direction [m]
// @param yDirLength  Total length in the y direction [m]
//**************************************************************************
modelMesh::modelMesh(int xCells, int yCells, double xDirLength, double
    yDirLength){
  numOfxCells = xCells;
  numOfyCells = yCells;
  xLength = xDirLength;
  yLength = yDirLength;
  dx = xLength/(float)numOfxCells;
  dy = yLength/(float)numOfyCells;
  numOfTotalCells = numOfxCells*numOfyCells;

  // Builds the geometry
  buildGeometry();
}
//*****************************************************************************
// Builds the problem mesh
//*****************************************************************************
void modelMesh::buildGeometry(){

  // Creates the mesh cells
  createCells();

  // Creates the cell faces
  createCellFaces();

  // Connects the mesh cells
  connectCells();
}

//*****************************************************************************
// Creates the nodes
//*****************************************************************************
void modelMesh::createCells(){
  int absIndex = 0;

  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      double x = i*dx + dx/2.; // x-coordinate
      double y = j*dy + dy/2.; // y-coordinate

      meshCell cell(i, j, absIndex, x, y,  dx, dy);
      meshCells.push_back(cell);
      absIndex++;
    }
  }
}
//*****************************************************************************
// Creates the cell faces
//*****************************************************************************
void modelMesh::createCellFaces(){
  int absIndex = 0, jmax;

  for (int i = 0; i <= 2*numOfxCells; i++){
    if (i%2){
      jmax = numOfyCells+1;
    }
    else{
      jmax = numOfyCells;
    }
    for (int j = 0; j < jmax; j++){
      meshCellFace face(i, j, absIndex);
      meshCellFaces.push_back(face);
      absIndex++;
    }
  }

}

//*****************************************************************************
// Connects the nodes
//*****************************************************************************
void modelMesh::connectCells(){
  meshCell* cellPtr;
  int kEast, kWest, kSouth, kNorth;
  double dxFlowArea, dyFlowArea;

  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      cellPtr = getCellByLoc(i,j);

      // Sets pointer to the cells
      meshCell* northCellPtr = getCellByLoc(i,j+1);
      meshCell* southCellPtr = getCellByLoc(i,j-1);
      meshCell* westCellPtr = getCellByLoc(i-1,j);
      meshCell* eastCellPtr = getCellByLoc(i+1,j);

      // Gets the absolute indices for the cell faces
      kEast = i*(2*numOfyCells+1) + j;
      kWest = (i+1)*(2*numOfyCells+1) + j;
      kSouth = numOfyCells*(2*i+1) + i + j;
      kNorth = numOfyCells*(2*i+1) + i + j + 1;

      // Sets the pointers for cell faces
      meshCellFace* eastFacePtr = &meshCellFaces[kEast];
      meshCellFace* westFacePtr = &meshCellFaces[kWest];
      meshCellFace* southFacePtr = &meshCellFaces[kSouth];
      meshCellFace* northFacePtr = &meshCellFaces[kNorth];

      // Sets the areas for flow
      dxFlowArea = (cellPtr->dx) ? cellPtr->dx : 1;
      dyFlowArea = (cellPtr->dy) ? cellPtr->dy : 1;

      connection northCon(northCellPtr, northFacePtr,
        dxFlowArea, -1., 0, cellPtr->dy);
      connection southCon(southCellPtr, southFacePtr,
        dxFlowArea, 1., 1, cellPtr->dy);
      connection eastCon(eastCellPtr, eastFacePtr,
        dyFlowArea, -1., 2, cellPtr->dx);
      connection westCon(westCellPtr, westFacePtr,
        dyFlowArea, 1., 3, cellPtr->dx);

      cellPtr->connections.push_back(northCon);
      cellPtr->connections.push_back(southCon);
      cellPtr->connections.push_back(eastCon);
      cellPtr->connections.push_back(westCon);
    }
  }
}

//*****************************************************************************
// Returns a pointer to the node by i,j location
//
// @param i  x index of cell
// @param j y index of cell
//*****************************************************************************
meshCell* modelMesh::getCellByLoc(int i, int j){
  meshCell *returnPtr;

  if (checkCellLoc(i, j)){
    int k = j + i*numOfyCells;
    assert(k <= meshCells.size() and k >= 0);
    returnPtr = &meshCells[k];
  }
  else{
    returnPtr = nullptr;
  }
  return returnPtr;
}

//*****************************************************************************
// Returns a pointer to the node by absolution index
//
// @param k  Absolute index of cell
//*****************************************************************************
meshCell* modelMesh::getCellByLoc(int k){
  meshCell *returnPtr;

  if (k <= meshCells.size() and k >= 0){
    returnPtr = &meshCells[k];
  }
  else{
    returnPtr = nullptr;
  }
  return returnPtr;
}

//*****************************************************************************
// Checks to see if i,j location is valid
//
// @param i  x index of cell
// @param j  y index of cell
//*****************************************************************************
bool modelMesh::checkCellLoc(int i, int j){

  if (i < 0)
    return false;
  else if (i >= numOfxCells)
    return false;
  else if (j < 0)
    return false;
  else if (j >= numOfyCells)
    return false;
  else
    return true;
}

//*****************************************************************************
// Sets a constant x velocity across the whole problem
//
// @param velocity  x velocity [m/s]
//*****************************************************************************
void modelMesh::setConstantXVelocity(double velocity){
  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      meshCell* cell = getCellByLoc(i,j);

      // loop over connections
      for (int conCount = 0; conCount < cell->connections.size(); conCount ++){
        connection thisCon = cell->connections[conCount];
        // sets the face velocity for east and west faces
        if (thisCon.loc == 2 or thisCon.loc == 3){
          thisCon.connectionFacePtr->vl = velocity;
        }
      }
    }
  }
}

//*****************************************************************************
// Sets a constant x velocity across a channel
//
// @param velocity  x velocity [m/s]
//*****************************************************************************
void modelMesh::setConstantXVelocity(double velocity, int row){
  for (int i = 0; i < numOfxCells; i++){
    meshCell* cell = getCellByLoc(i,row);

    // loop over connections
    for (int conCount = 0; conCount < cell->connections.size(); conCount ++){
      connection thisCon = cell->connections[conCount];
      // sets the face velocity for east and west faces
      if (thisCon.loc == 2 or thisCon.loc == 3){
        thisCon.connectionFacePtr->vl = velocity;
      }
    }
  }
}

//*****************************************************************************
// Sets a constant y velocity across the whole problem
//
// @param velocity  y velocity [m/s]
//*****************************************************************************
void modelMesh::setConstantYVelocity(double velocity){
  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      meshCell* cell = getCellByLoc(i,j);

      // loop over connections
      for (int conCount = 0; conCount < cell->connections.size(); conCount ++){
        connection thisCon = cell->connections[conCount];
        // sets the face velocity for north and south faces
        if (thisCon.loc == 0 or thisCon.loc == 1){
          thisCon.connectionFacePtr->vl = velocity;
        }
      }
    }
  }
}

//*****************************************************************************
// Sets a constant y velocity across a channel
//
// @param velocity  y velocity [m/s]
//*****************************************************************************
void modelMesh::setConstantYVelocity(double velocity, int column){
  for (int j = 0; j < numOfyCells; j++){
    meshCell* cell = getCellByLoc(column,j);

    // loop over connections
    for (int conCount = 0; conCount < cell->connections.size(); conCount ++){
      connection thisCon = cell->connections[conCount];
      // sets the face velocity for north and south faces
      if (thisCon.loc == 0 or thisCon.loc == 1){
        thisCon.connectionFacePtr->vl = velocity;
      }
    }
  }
}

//*****************************************************************************
// Sets the temperature in the whole problem
//
// @param temp    Temperature in kelvin
//*****************************************************************************
void modelMesh::setSystemTemperature(double temp){
  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      meshCell* cell = getCellByLoc(i,j);
      cell->setTemperature(temp);
    }
  }
}

//*****************************************************************************
// Sets the temperature in the whole problem from a csv file
// The order of the values in the csv file is very important. For a 3x4 example
// with 3 cells in the x-direction and 4 in the y-direction the format should
// be
//
//    0   1   2
// 0 v1  v2  v3
// 1 v4  v5  v6
// 2 v7  v8  v9
// 3 v10 v11 v12
//
// Temperature in kelvin
//
// @param path  Path to the csv file
//*****************************************************************************
void modelMesh::setSystemTemperature(const std::string path){
  // Init the matrix to hold the temp values
  MatrixD tempatureMat = MatrixD::Zero(numOfyCells, numOfxCells);
  // read in the csv file
  readCSV(tempatureMat, path);

  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      double temp = tempatureMat(j,i);
      meshCell* cell = getCellByLoc(i,j);
      cell->setTemperature(temp);
    }
  }
}

//*****************************************************************************
// Sets the pressure in the whole problem
//
// @param pressure    Pressure in Pa
//*****************************************************************************
void modelMesh::setSystemPressure(double pressure){
  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      meshCell* cell = getCellByLoc(i,j);
      cell->setPressure(pressure);
    }
  }
}

//*****************************************************************************
// Sets the pressure in the whole problem from a csv file
// The order of the values in the csv file is very important. For a 3x4 example
// with 3 cells in the x-direction and 4 in the y-direction the format should
// be
//
//    0   1   2
// 0 v1  v2  v3
// 1 v4  v5  v6
// 2 v7  v8  v9
// 3 v10 v11 v12
//
// Presure in Pa
//
// @param path  Path to the csv file
//*****************************************************************************
void modelMesh::setSystemPressure(const std::string path){
  // Init the matrix to hold the temp values
  MatrixD pressureMat = MatrixD::Zero(numOfyCells, numOfxCells);
  // read in the csv file
  readCSV(pressureMat, path);

  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      double pres = pressureMat(j,i);
      meshCell* cell = getCellByLoc(i,j);
      cell->setPressure(pres);
    }
  }
}

//*****************************************************************************
// Sets the neutron flux in the whole problem
//
// @param phi    Neutron flux in 1/cm^2/s
//*****************************************************************************
void modelMesh::setSystemNeutronFlux(double phi){
  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      meshCell* cell = getCellByLoc(i,j);
      cell->setNeutronFlux(phi);
    }
  }
}

//*****************************************************************************
// Sets the neutron flux in the whole problem from a csv file
// The order of the values in the csv file is very important. For a 3x4 example
// with 3 cells in the x-direction and 4 in the y-direction the format should
// be
//
//    0   1   2
// 0 v1  v2  v3
// 1 v4  v5  v6
// 2 v7  v8  v9
// 3 v10 v11 v12
//
// neutron flux in 1/cm^2/s
//
// @param path  Path to the csv file
//*****************************************************************************
void modelMesh::setSystemNeutronFlux(const std::string path){
  // Init the matrix to hold the temp values
  MatrixD neutronMat = MatrixD::Zero(numOfyCells, numOfxCells);
  // read in the csv file
  readCSV(neutronMat, path);

  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      double phi = neutronMat(j,i);
      meshCell* cell = getCellByLoc(i,j);
      cell->setNeutronFlux(phi);
    }
  }
}

//*****************************************************************************
// Sets the gas phase interfacial area concentration
//
// @param intAreaCon  Interfacial area concentration [1/m]
//*****************************************************************************
void modelMesh::setSystemGasInterfacialAreaCon(double intAreaCon){
  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      meshCell* cell = getCellByLoc(i,j);
      cell->setGasInterfacialAreaCon(intAreaCon);
    }
  }
}

//*****************************************************************************
// Sets the gas phase interfacial area concentration in the whole problem
// from a csv file
// The order of the values in the csv file is very important. For a 3x4 example
// with 3 cells in the x-direction and 4 in the y-direction the format should
// be
//
//    0   1   2
// 0 v1  v2  v3
// 1 v4  v5  v6
// 2 v7  v8  v9
// 3 v10 v11 v12
//
// intAreaCon in 1/m
//
// @param path  Path to the csv file
//*****************************************************************************
void modelMesh::setSystemGasInterfacialAreaCon(const std::string path){
  // Init the matrix to hold the values
  MatrixD gasIntAreaMat = MatrixD::Zero(numOfyCells, numOfxCells);
  // read in the csv file
  readCSV(gasIntAreaMat, path);

  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      double intAreaCon = gasIntAreaMat(j,i);
      meshCell* cell = getCellByLoc(i,j);
      cell->setGasInterfacialAreaCon(intAreaCon);
    }
  }
}

//*****************************************************************************
// Sets the surface interfacial area concentration
//
// @param intAreaCon  Interfacial area concentration [1/m]
//*****************************************************************************
void modelMesh::setSystemWallInterfacialAreaCon(double intAreaCon){
  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      meshCell* cell = getCellByLoc(i,j);
      cell->setWallInterfacialAreaCon(intAreaCon);
    }
  }
}

//*****************************************************************************
// Sets the wall interfacial area concentration in the whole problem
// from a csv file
// The order of the values in the csv file is very important. For a 3x4 example
// with 3 cells in the x-direction and 4 in the y-direction the format should
// be
//
//    0   1   2
// 0 v1  v2  v3
// 1 v4  v5  v6
// 2 v7  v8  v9
// 3 v10 v11 v12
//
// intAreaCon in 1/m
//
// @param path  Path to the csv file
//*****************************************************************************
void modelMesh::setSystemWallInterfacialAreaCon(const std::string path){
  // Init the matrix to hold the values
  MatrixD wallIntAreaMat = MatrixD::Zero(numOfyCells, numOfxCells);
  // read in the csv file
  readCSV(wallIntAreaMat, path);

  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      double wallAreaCon = wallIntAreaMat(j,i);
      meshCell* cell = getCellByLoc(i,j);
      cell->setWallInterfacialAreaCon(wallAreaCon);
    }
  }
}

//*****************************************************************************
// Sets the gas void fraction
//
// @param fract  Gas void fraction
//*****************************************************************************
void modelMesh::setSystemGasVoidFraction(double fract){
  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      meshCell* cell = getCellByLoc(i,j);
      cell->setGasVoidFraction(fract);
    }
  }
}

//*****************************************************************************
// Sets the void fraction in the whole problem
// from a csv file
// The order of the values in the csv file is very important. For a 3x4 example
// with 3 cells in the x-direction and 4 in the y-direction the format should
// be
//
//    0   1   2
// 3 v1  v2  v3
// 2 v4  v5  v6
// 1 v7  v8  v9
// 0 v10 v11 v12
//
// fract in void fraction
//
// @param path  Path to the csv file
//*****************************************************************************
void modelMesh::setSystemGasVoidFraction(const std::string path){
  // Init the matrix to hold the values
  MatrixD voidMat = MatrixD::Zero(numOfyCells, numOfxCells);
  // read in the csv file
  readCSV(voidMat, path);

  for (int i = 0; i < numOfxCells; i++){
    for (int j = 0; j < numOfyCells; j++){
      double fract = voidMat(j,i);
      meshCell* cell = getCellByLoc(i,j);
      cell->setGasVoidFraction(fract);
    }
  }
}
//*****************************************************************************
// Set system variable value
//
// @param variable  Name of the variable to set
// @param value     Variable value
//*****************************************************************************
void modelMesh::setSystemVariableValue(const std::string variable, const double value){
  if (variable == "x_velocity"){
    setConstantXVelocity(value);
  }
  else if (variable == "y_velocity"){
    setConstantYVelocity(value);
  }
  else if (variable == "temperature"){
    setSystemTemperature(value);
  }
  else if (variable == "pressure"){
    setSystemPressure(value);
  }
  else if (variable == "neturon_flux"){
    setSystemNeutronFlux(value);
  }
  else if (variable == "gas_interfacial_area_concentration"){
    setSystemGasInterfacialAreaCon(value);
  }
  else if (variable == "wall_interfacial_area_concentration"){
    setSystemWallInterfacialAreaCon(value);
  }
  else if (variable == "gas_volume_fraction"){
    setSystemGasVoidFraction(value);
  }
  else{
    std::string errorMessage = "Variable " + variable + " is not a supported type.";
    libowskiException::runtimeError(errorMessage);
  }
}

//*****************************************************************************
// Set cell temperature
//
// @param i      x cell index
// @param j      y cell index
// @param temp    Temperature in kelvin
//*****************************************************************************
void modelMesh::setCellTemperature(int i, int j, double temp){
  meshCell* cell = getCellByLoc(i,j);
  cell->setTemperature(temp);
}

//*****************************************************************************
// Set cell pressure
//
// @param i        x cell index
// @param j        y cell index
// @param Pressure  Pressure in Pa
//*****************************************************************************
void modelMesh::setCellPressure(int i, int j, double pressure){
  meshCell* cell = getCellByLoc(i,j);
  cell->setPressure(pressure);
}

//*****************************************************************************
// Set cell neutron flux
//
// @param i        x cell index
// @param j        y cell index
// @param phi      Neutron flux in 1/cm^2/s
//*****************************************************************************
void modelMesh::setCellNeutronFlux(int i, int j, double phi){
  meshCell* cell = getCellByLoc(i,j);
  cell->setNeutronFlux(phi);
}

//*****************************************************************************
// Set cell interfacial area conentration for gas phase
//
// @param i        x cell index
// @param j        y cell index
// @param intAreaCon  Interfacial area concentration 1/m
//*****************************************************************************
void modelMesh::setCellGasInterfacialAreaCon(int i, int j, double intAreaCon){
  meshCell* cell = getCellByLoc(i,j);
  cell->setGasInterfacialAreaCon(intAreaCon);
}

//*****************************************************************************
// Set cell surface interfacial area conentration
//
// @param i        x cell index
// @param j        y cell index
// @param intAreaCon  Interfacial area concentration 1/m
//*****************************************************************************
void modelMesh::setCellWallInterfacialAreaCon(int i, int j, double intAreaCon){
  meshCell* cell = getCellByLoc(i,j);
  cell->setWallInterfacialAreaCon(intAreaCon);
}

//*****************************************************************************
// Set cell gas void fraction
//
// @param i      x cell index
// @param j      y cell index
// @param fract  Gas void fraction
//*****************************************************************************
void modelMesh::setGasVoidFraction(int i, int j, double fract){
  meshCell* cell = getCellByLoc(i,j);
  cell->setGasVoidFraction(fract);
}

//*****************************************************************************
// Adds a surface to a cell
//
// @param i        x cell index
// @param j        y cell index
// @param loc      Location of the surface in the cell
//              north
//              south
//              east
//              west
//*****************************************************************************
void modelMesh::addSurface(int i, int j, std::string loc){
  int locID = -1;
  meshCell* cell = getCellByLoc(i,j);

  if (loc == "north"){locID = 0;};
  if (loc == "south"){locID = 1;};
  if (loc == "east"){locID = 2;};
  if (loc == "west"){locID = 3;};
  assert(locID != -1);

  cell->addSurface(locID);
}

//*****************************************************************************
// Adds a surface to the mesh boundary
//
// @param loc      Location of the surface boundary
//              north
//              south
//              east
//              west
//*****************************************************************************
void modelMesh::addBoundarySurface(std::string loc){
  int xCellMax = numOfxCells - 1;
  int xCellMin = 0;
  int yCellMax = numOfyCells - 1;
  int yCellMin = 0;
  int locID = -1;

  if (loc == "north"){locID = 0;};
  if (loc == "south"){locID = 1;};
  if (loc == "east"){locID = 2;};
  if (loc == "west"){locID = 3;};
  assert(locID != -1);

  switch(locID){

    // North location
    case 0: {
      for (int i = 0; i < numOfxCells; i++){
        addSurface(i, yCellMax, "north");
      }
      break;
    }
    // South location
    case 1: {
      for (int i = 0; i < numOfxCells; i++){
        addSurface(i, yCellMin, "south");
      }
      break;
    }
    // East location
    case 2: {
      for (int j = 0; j < numOfyCells; j++){
        addSurface(xCellMax, j, "east");
      }
      break;
    }
    // West location
    case 3: {
      for (int j = 0; j < numOfyCells; j++){
        addSurface(xCellMin, j, "west");
      }
      break;
    }
  }

}

//*****************************************************************************
// Cleans the model
//*****************************************************************************
void modelMesh::clean(){

  xLength = 0.0;
  yLength = 0.0;
  numOfxCells = 0;
  numOfyCells = 0;
  numOfTotalCells = 0;
  meshCells.clear();
  meshCellFaces.clear();
  dx = 0.0;
  dy = 0.0;
}
