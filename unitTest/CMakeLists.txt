project(unitTest)

# The test program
add_executable(unitTestSolver testSolver.cpp)
add_executable(unitTestMesh testMesh.cpp)
if(${MPI_FOUND})
	# Solver test
	target_link_libraries(unitTestSolver PRIVATE liblibowski ${MPI_LIBRARIES})
	target_compile_definitions(unitTestSolver PRIVATE HAVE_MPI)
	add_test(test_solver ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 unitTestSolver)

	# Model mesh test
	target_link_libraries(unitTestMesh PRIVATE liblibowski ${MPI_LIBRARIES})
	target_compile_definitions(unitTestMesh PRIVATE HAVE_MPI)
	add_test(test_mesh ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 unitTestMesh)
else()
	# Solver test
	target_link_libraries(unitTestSolver PRIVATE liblibowski)
	add_test(test_Solver unitTestSolver)

	# Model mesh test
	target_link_libraries(unitTestMesh PRIVATE liblibowski)
	add_test(test_mesh unitTestMesh)
endif()
