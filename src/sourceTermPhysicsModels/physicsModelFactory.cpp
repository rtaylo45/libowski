#include "physicsModelFactory.h"

//**************************************************************************
// Methods for the factory
//**************************************************************************
physicsModel *physicsModelFactory::getPhysicsModel(std::string type){
	physicsModel *physicModel = nullptr;

	if (type == "generic"){
		physicModel = new generic();
		return physicModel;
	}

	return physicModel;

}
