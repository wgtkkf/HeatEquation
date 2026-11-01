#
cmake_minimum_required(VERSION 3.20)

#
project(ConductionFEM CXX)

#
set(CMAKE_CXX_STANDARD 20) # C++20
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 1. Bundle all core logic into a static library called 'FEMCore'
add_library(FEMCore STATIC
    src/Parameter.cpp
    src/Mesh.cpp
    src/Boundary.cpp
)

# Make sure the library knows where the header files are
target_include_directories(FEMCore PUBLIC include)

# 2. My main executable
add_executable(bmt src/source.cpp)

# Link the core library to my executable
target_link_libraries(bmt PRIVATE FEMCore)

# ---------------------------------------------------
# TESTING SECTION
# ---------------------------------------------------

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
add_executable(run_tests
    tests/test_mesh.cpp
)

# Link my test executable to the Google test libraries
target_link_libraries(run_tests PRIVATE FEMCore GTest::gtest_main)

# Register the test with CMake
add_test(NAME UnitTests COMMAND run_tests)

# ---------------------------------------------------
# Memo
# build commands:
# 0. cd build
# 1. cmake -DCMAKE_CXX_COMPILER=g++-13 ..    
# 2. make
# 3. ./bmt
# ---------------------------------------------------
