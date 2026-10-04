#!/bin/bash

echo "Cleaning old build..."
rm -rf build/

echo "Configuring with GCC-13..."
cmake -B build -DCMAKE_CXX_COMPILER=g++-13

echo "Compiling..."
cmake --build build
