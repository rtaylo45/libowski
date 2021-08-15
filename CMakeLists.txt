cmake_minimum_required(VERSION 3.13.0)
project(libowski)

# Sets directories
set(SRC_DIR src)
set(LIB_DIR lib)
set(UTIL_DIR ${LIB_DIR}/utils)
set(MASSTRANSFER_DIR ${SRC_DIR}/massTransferModels)

set(EIGEN3_INCLUDE_DIR ${LIB_DIR}/eigen3)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(MPI)
if(${MPI_FOUND})
	include_directories(${MPI_INCLUDE_PATH})
endif()

if(${MPI_FOUND})
	set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
	set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
endif()

# Links the file system library 
link_libraries(stdc++fs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3")

# util directory
add_subdirectory(${UTIL_DIR})
# Add the source directory
add_subdirectory(${SRC_DIR})
# Add unit test directory
add_subdirectory(unitTest)
# Add case studie directory
add_subdirectory(caseStudies)

# Enable testing
ENABLE_TESTING()
