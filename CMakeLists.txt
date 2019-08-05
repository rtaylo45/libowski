cmake_minimum_required(VERSION 3.0)
project(libowski)

# Sets source directory
set(SRC_DIR src)
set(LIB_DIR lib)


# The find cmake flag. I have included eigen3 with my softwar now
# so the user doesn't need to install it.
#find_package(Eigen3 REQUIRED)
# Path to eigen dir
#set(${EIGEN3_INCLUDE_DIR})

set(EIGEN3_INCLUDE_DIR ${LIB_DIR}/eigen3)
include_directories(${EIGEN3_INCLUDE_DIR})

# Sets the standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")

# Project sources
set(PROJECT_SOURCES ${SRC_DIR}/coreSolver.h
					${SRC_DIR}/coreSolver.cpp
)

# main excutable
add_executable(
	libowski
	${SRC_DIR}/main.cpp
	${PROJECT_SOURCES}
)

# unit test
add_executable(
	test
	${SRC_DIR}/main.cpp
	${PROJECT_SOURCES}
)


