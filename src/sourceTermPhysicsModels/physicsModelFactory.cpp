#include "physicsModelFactory.h"

//**************************************************************************
// Methods for the factory
//**************************************************************************
physicsModels *physicsModelFactory::getPhysicsModel(std::string type){
	physicsModels *physicModel = nullptr;

	if (type == "generic"){
		//physicModel = new generic(meshCell* species*);
		return physicModel;
	}

	return physicModel;

}
