project(caseStudies)

add_executable(caseStudies caseStudies.cpp)

if(${MPI_FOUND})
	target_link_libraries(caseStudies PRIVATE liblibowski ${MPI_LIBRARIES})
	target_compile_definitions(caseStudies PRIVATE HAVE_MPI)
else()
	target_link_libraries(caseStudies PRIVATE liblibowski)
endif()
