project(unitTest)

# The test program
add_executable(unitTest test.cpp)
if(${MPI_FOUND})
	target_link_libraries(unitTest PRIVATE liblibowski ${MPI_LIBRARIES})
	target_compile_definitions(unitTest PRIVATE HAVE_MPI)
	add_test(test_main ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 unitTest)
else()
	target_link_libraries(unitTest PRIVATE liblibowski)
	add_test(test_main unitTest)
endif()
