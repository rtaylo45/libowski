project(utils)

set(UTIL_SOURCES
  src/linearAlgebra.cpp
  src/exception.cpp
  src/utilBase.cpp
  src/sys.cpp
  src/Faddeeva.cc)

add_library(utillibowski STATIC ${UTIL_SOURCES})
#target_include_directories(utillibowski
#PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
