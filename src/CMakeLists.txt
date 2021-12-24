project(libowski_prog)

set(PROJECT_SOURCES
  mpiProcess.cpp
  modelMesh.cpp
  meshCell.cpp
  meshCellFace.cpp
  species.cpp
  speciesDriver.cpp
  convectionLimiter.cpp
  matrixExponential.cpp
  ODEintegrator.cpp
  cellConnection.cpp
  surface.cpp
  scalarData.cpp
)

# All sources that also need to be tested in unit tests go into a static library
add_library(liblibowski STATIC ${PROJECT_SOURCES})
#target_include_directories(liblibowski PUBLIC ${INCLUDE_DIR})
target_link_libraries(liblibowski PUBLIC utillibowski)

# Adds sub source directories
add_subdirectory(sourceTerms)

# The main program
add_executable(libowski main.cpp)
if(${MPI_FOUND})
  target_link_libraries(liblibowski PRIVATE ${MPI_LIBRARIES})
  target_link_libraries(libowski PRIVATE liblibowski ${MPI_LIBRARIES})
  target_compile_definitions(libowski PRIVATE HAVE_MPI)
  target_compile_definitions(liblibowski PRIVATE HAVE_MPI)
else()
  target_link_libraries(libowski PRIVATE liblibowski)
endif()
