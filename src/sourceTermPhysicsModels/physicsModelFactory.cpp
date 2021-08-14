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
	else if (type == "genericRemoval"){
		physicModel = new genericRemoval();
		return physicModel;
	}
	else if (type == "neutronInduced"){
		physicModel = new neutronInducedReactions();
		return physicModel;
	}
	else if (type == "wallDeposition"){
		physicModel = new wallDeposition();
		return physicModel;
	}
	else if (type == "gasSparging"){
		physicModel = new gasSparging();
		return physicModel;
	}
	else {
		std::string errorMessage =
			" You have selected a physics type that is not\n"
			" avaliable. Avaliable models are\n\n"
			" generic\n"
			" genericRemoval\n"
			" neutronInduced\n"
			" wallDeposition\n";
		libowskiException::runtimeError(errorMessage);
		return physicModel;
	}

	return physicModel;

}
