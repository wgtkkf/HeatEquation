# HeatEquation
Two-dimensional finite element modeling of heat equasion.

# A prerequisit for C-language code
Create two folders: vtk and excel

# Run
1. gcc source.c -o source -lm
2. ./source

# C++ code

```text
my_multiphysics_solver/
├── CMakeLists.txt       # The master instruction manual for how to build the code
├── Dockerfile           # The container configuration
├── include/             # ONLY header files (.h or .hpp) go here
│   ├── Mesh.hpp
│   └── HeatSolver.hpp
├── src/                 # ONLY implementation files (.cpp) go here
│   ├── main.cpp         # The entry point that runs the simulation
│   ├── Mesh.cpp
│   └── HeatSolver.cpp
├── tests/               # Your GoogleTest verification scripts go here
│   ├── CMakeLists.txt   # A sub-instruction file just for building tests
│   └── test_heat.cpp
└── build/               # The "trash can" folder