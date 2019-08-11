project(libowski_prog)

set(PROJECT_SOURCES coreSolver.h coreSolver.cpp)

# All sources that also need to be tested in unit tests go into a static library
add_library(liblibowski STATIC ${PROJECT_SOURCES})
message(${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(liblibowski PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# The main program
add_executable(libowski main.cpp)
target_link_libraries(libowski PRIVATE liblibowski ${MPI_LIBRARIES})
