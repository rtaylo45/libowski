project(caseStudies)

add_executable(caseStudie caseStudie.cpp)

if(${MPI_FOUND})
	target_link_libraries(caseStudie PRIVATE liblibowski ${MPI_LIBRARIES})
	target_compile_definitions(caseStudie PRIVATE HAVE_MPI)
else()
	target_link_libraries(caseStudie PRIVATE liblibowski)
endif()
