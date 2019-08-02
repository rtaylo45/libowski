cmake_minimum_required(VERSION 3.0)
project(libowski)

# Sets source directory
set(SRC_DIR src)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Sets the standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")

# adds the source files
add_executable(
	libowski
	${SRC_DIR}/main.cpp
	${SRC_DIR}/coreSolver.h
	${SRC_DIR}/coreSolver.cpp
)

