# CMake generated Testfile for 
# Source directory: /home/tubuntu/FEM_MyTrials/HeatEquation/cpp
# Build directory: /home/tubuntu/FEM_MyTrials/HeatEquation/cpp/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[UnitTests]=] "/home/tubuntu/FEM_MyTrials/HeatEquation/cpp/build/run_tests")
set_tests_properties([=[UnitTests]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/tubuntu/FEM_MyTrials/HeatEquation/cpp/CMakeLists.txt;49;add_test;/home/tubuntu/FEM_MyTrials/HeatEquation/cpp/CMakeLists.txt;0;")
subdirs("_deps/googletest-build")
