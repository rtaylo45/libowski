//*****************************************************************************
// Author: Zack Taylor
//
// Exception handeling 
//*****************************************************************************
#include "exception.h"

//*****************************************************************************
// Prints out runtime error header 
//*****************************************************************************
void libowskiException::runtimeErrorOutputHeader(){
	std::cout << "**************************************************\n"
					 " You have encountered a runtime error in libowski.\n"
					 " Your error message is:\n" << std::endl;
}

//*****************************************************************************
// Prints out runtime error butt
//*****************************************************************************
void libowskiException::runtimeErrorOutputButt(){
	std::cout <<"**************************************************"<< std::endl;
}
//*****************************************************************************
// Throw runtime error and kill program
//
// @param errorMessage	The error message shown to the terminal screen
//*****************************************************************************
void libowskiException::runtimeError(const std::string errorMessage){

	runtimeErrorOutputHeader();
	std::cout<<errorMessage << std::endl;
	runtimeErrorOutputButt();
	exit(1);
}
