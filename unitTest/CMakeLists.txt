project(unitTest)

# The test program
add_executable(unitTestMatExpSolver testMatExpSolver.cpp)
add_executable(unitTestIntegratorSolvers testIntegratorSolvers.cpp)
add_executable(unitTestMesh testMesh.cpp)
add_executable(unitTestVerification verificationProblems.cpp)
add_executable(unitTestUtils testUtils.cpp)
add_executable(unitTestMoleProblems testMoleProblems.cpp)
add_executable(unitTestSpatialConvergence testSpatialConvergence.cpp)

# Unit test path
set(TEST_PATH ${CMAKE_CURRENT_BINARY_DIR})
if(${PARALLEL})
  # Integrator Solver test
  target_link_libraries(unitTestIntegratorSolvers PRIVATE liblibowski ${MPI_LIBRARIES})
  target_compile_definitions(unitTestIntegratorSolvers PRIVATE HAVE_MPI)
  add_test(test_Integrator_solver ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${TEST_PATH}/unitTestIntegratorSolvers)

  # Matrix Exponential Solver test
  target_link_libraries(unitTestMatExpSolver PRIVATE liblibowski ${MPI_LIBRARIES})
  target_compile_definitions(unitTestMatExpSolver PRIVATE HAVE_MPI)
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

  # mole problems
  target_link_libraries(unitTestMoleProblems PRIVATE liblibowski ${MPI_LIBRARIES})
  target_compile_definitions(unitTestMoleProblems PRIVATE HAVE_MPI)
  add_test(test_Mole ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${TEST_PATH}/unitTestMoleProblems)

  # Spatial convergence problems
  target_link_libraries(unitTestSpatialConvergence PRIVATE liblibowski ${MPI_LIBRARIES})
  target_compile_definitions(unitTestSpatialConvergence PRIVATE HAVE_MPI)
  add_test(test_Spatial_convergence ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${TEST_PATH}/unitTestSpatialConvergence)

else()
  # Integrator Solver test
  target_link_libraries(unitTestIntegratorSolvers PRIVATE liblibowski)
  add_test(test_Integrator_solvers ${TEST_PATH}/unitTestIntegratorSolvers)

  # Matrix exp Solver test
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

  # Mole problems
  target_link_libraries(unitTestMoleProblems PRIVATE liblibowski)
  add_test(test_Mole ${TEST_PATH}/unitTestMoleProblems)

  # Flux convergence problems
  target_link_libraries(unitTestSpatialConvergence PRIVATE liblibowski)
  add_test(test_Spatial_convergence ${TEST_PATH}/unitTestSpatialConvergence)
endif()
