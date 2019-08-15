project(unitTest)

# The test program
add_executable(unitTest test.cpp)
target_link_libraries(unitTest PRIVATE liblibowski ${MPI_LIBRARIES})
add_test(test_main ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 unitTest)
