//*****************************************************************************
// Author: Zack Taylor
//
// Exception handeling class for libowski
//*****************************************************************************
#ifndef LIBOWSKIEXCEPTION_H
#define LIBOWSKIEXCEPTION_H
#include <stdexcept>
#include <string>
#include <iostream>

//*****************************************************************************
// Exception handeling class
//*****************************************************************************
class libowskiException{
	private:
	//**************************************************************************
	// Error header for runtime error
	//**************************************************************************
	static void runtimeErrorOutputHeader();
	//**************************************************************************
	// Error butt for runtime error
	//**************************************************************************
	static void runtimeErrorOutputButt();
	public:
	//**************************************************************************
	// Runtime error 
	//**************************************************************************
	static void runtimeError(const std::string&);
};

#endif
