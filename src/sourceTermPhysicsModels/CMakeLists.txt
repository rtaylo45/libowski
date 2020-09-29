target_sources(liblibowski
	PUBLIC
	physicsModelABC.cpp
	physicsModelFactory.cpp
	generic.cpp
)
target_include_directories(liblibowski PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
