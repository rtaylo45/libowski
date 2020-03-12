project(unitTest)

# The test program
add_executable(unitTestMatExpSolver testMatExpSolver.cpp)
add_executable(unitTestMesh testMesh.cpp)
add_executable(unitTestVerification verificationProblems.cpp)
add_executable(unitTestUtils testUtils.cpp)

# Unit test path
set(TEST_PATH ${CMAKE_CURRENT_BINARY_DIR})
if(${MPI_FOUND})
	# Matrix Exponential Solver test
	target_link_libraries(unitTestMatExpSolver PRIVATE liblibowski ${MPI_LIBRARIES})
	target_compile_definitions(unitTestSolver PRIVATE HAVE_MPI)
	add_test(test_MatExp_solver ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${TEST_PATH}/unitTestMatExpSolver)

	# Model mesh test
	target_link_libraries(unitTestMesh PRIVATE liblibowski ${MPI_LIBRARIES})
	target_compile_definitions(unitTestMesh PRIVATE HAVE_MPI)
	add_test(test_mesh ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${TEST_PATH}/unitTestMesh)

	# Verification problem
	target_link_libraries(unitTestVerification PRIVATE liblibowski ${MPI_LIBRARIES})
	target_compile_definitions(unitTestVerification PRIVATE HAVE_MPI)
	add_test(test_Verification ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${TEST_PATH}/unitTestVerification)

	# Test Utilities
	target_link_libraries(unitTestUtils PRIVATE liblibowski ${MPI_LIBRARIES})
	target_compile_definitions(unitTestUtils PRIVATE HAVE_MPI)
	add_test(test_Utils ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${TEST_PATH}/unitTestUtils)
else()
	# Solver test
	target_link_libraries(unitTestMatExpSolver PRIVATE liblibowski)
	add_test(test_MatExp_solver ${TEST_PATH}/unitTestMatExpSolver)

	# Model mesh test
	target_link_libraries(unitTestMesh PRIVATE liblibowski)
	add_test(test_mesh ${TEST_PATH}/unitTestMesh)

	# Verification problem
	target_link_libraries(unitTestVerification PRIVATE liblibowski)
	add_test(test_Verification ${TEST_PATH}/unitTestVerification)

	# Test Utilities
	target_link_libraries(unitTestUtils PRIVATE liblibowski)
	add_test(test_Utils ${TEST_PATH}/unitTestUtils)
endif()
