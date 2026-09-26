#
cmake_minimum_required(VERSION 3.0)

#
project(ConductionFEM CXX)

#
set(CMAKE_CXX_STANDARD 20) # C++20
set(CMAKE_CXX_STANDARD_REQUIRED True)

#
add_executable(bmt
    src/source.cpp
    src/Parameter.cpp
)

#
target_include_directories(bmt PRIVATE include)

# memo
# build commands:
# 0. cd build
# 1. cmake -DCMAKE_CXX_COMPILER=g++-13 ..    
# 2. make
# 3. ./bmt
