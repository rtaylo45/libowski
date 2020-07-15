//*****************************************************************************
// Author: Zack Taylor
//
// A collection of functions used for system information/stuff
//*****************************************************************************
#include <fstream>
#include <string>
#include <filesystem>
#include <exception.h>

//*****************************************************************************
// Gets the path for the data folder of libowski
//*****************************************************************************
std::string getDataPath();

//*****************************************************************************
// Checks to see if a file exist or not
//*****************************************************************************
void checkFileExists(std::string fname);
