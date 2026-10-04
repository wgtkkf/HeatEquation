#
cmake_minimum_required(VERSION 3.20)

#
project(ConductionFEM CXX)

#
set(CMAKE_CXX_STANDARD 20) # C++20
set(CMAKE_CXX_STANDARD_REQUIRED True)

#
add_executable(bmt
    src/source.cpp
    src/Parameter.cpp
    src/Mesh.cpp
)

#
target_include_directories(bmt PRIVATE include)

#########################################
# memo
# build commands:
# 0. cd build
# 1. cmake -DCMAKE_CXX_COMPILER=g++-13 ..    
# 2. make
# 3. ./bmt
#########################################

## enable CTest
enable_testing()

# Fetch Google Test automatically
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)

# Create the test executable by compiling Mesh.cpp with my test file
add_executable(run_tests src/Mesh.cpp tests/test_mesh.cpp)

# Link my test executable to the Google test libraries
target_link_libraries(run_tests GTest::gtest_main)

# Register the test with CMake
add_test(NAME UnitTests COMMAND run_tests)
