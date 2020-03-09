project(utils)

set(UTIL_SOURCES linearAlgebra.cpp exception.cpp utilBase.cpp)

add_library(utillibowski STATIC ${UTIL_SOURCES})
target_include_directories(utillibowski PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
message(${CMAKE_CURRENT_SOURCE_DIR})
