project(unitTest)

# The test program
add_executable(unitTest test.cpp)
target_link_libraries(unitTest PRIVATE liblibowski)
add_test(test_main unitTest)
