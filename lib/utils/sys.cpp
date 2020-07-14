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
