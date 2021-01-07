project(utils)

set(UTIL_SOURCES linearAlgebra.cpp exception.cpp utilBase.cpp sys.cpp Faddeeva.cc)

add_library(utillibowski STATIC ${UTIL_SOURCES})
target_include_directories(utillibowski PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
