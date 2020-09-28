//*****************************************************************************
// Author: Zack Taylor
//
//*****************************************************************************
#ifndef PHYSICSMODELFACTORY_H
#define PHYSICSMODELFACTORY_H
#include <string>
#include "physicsModelABC.h"

//*****************************************************************************
// Physics model factory
//*****************************************************************************
class	physicsModelFactory{
	public:
	static physicsModels * getPhysicsModel(std::string);

};

#endif
