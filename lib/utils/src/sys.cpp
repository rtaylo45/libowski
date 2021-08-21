//*****************************************************************************
// Author: Zack Taylor
//
// A collection of functions used for system information/stuff
//*****************************************************************************
#include "sys.h"

namespace fs = std::filesystem;

//*****************************************************************************
// Gets the path for the data folder of libowski
//*****************************************************************************
std::string getDataPath(){
   std::string dataPath, estr;
   fs::path p = fs::current_path();
   for(auto& e : p){
      dataPath += e;
      if (e == "libowski"){
         break;
      }
      if (e != "/"){
         dataPath += "/";
      }
   }
   return dataPath +"/src/data/";

}

//*****************************************************************************
// Checks to see if a file exist or not
//
// @param fname	absolute location of the file
//*****************************************************************************
void checkFileExists(std::string fname){
	if (not fs::exists(fname)){
	   std::string errorMessage =
	      " Unable to find the file " + fname +"\n";
	   libowskiException::runtimeError(errorMessage);
	}
}
