//*****************************************************************************
// Author: Zack Taylor
//
//*****************************************************************************
#ifndef PHYSICSMODELFACTORY_H
#define PHYSICSMODELFACTORY_H
#include <string>
#include "physicsModelABC.h"
#include "generic.h"

//*****************************************************************************
// Physics model factory
//*****************************************************************************
class	physicsModelFactory{
	public:
	static physicsModel *getPhysicsModel(std::string);

};

#endif
