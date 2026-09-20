## About this repository
Two-dimensional finite element modeling of heat equasion.

## Boundary conditions
Fixed temperature, a schematic is in preparation.

## Analytical solution
An analytical solution was given by a python script.

## A prerequisit for C-language code
Please create two folders: vtk and excel

## Run C-language code
1 gcc source.c -o source -lm  
2 ./source

## C++ version
```text
cpp_version/
├── CMakeLists.txt       # 
├── Dockerfile           # 
├── include/             # .hpp files
│   ├── Timer.hpp
│   ├── Inputs.hpp
│   ├── Mesh.hpp
│   ├── Boundary.hpp
│   ├── Matrix.hpp
│   ├── Solver.hpp      
│   └── Outputs.hpp      
├── src/                 # .cpp files
│   ├── main.cpp         # entry point
│   ├── Timer.cpp
│   ├── Inputs.cpp
│   ├── Mesh.cpp
│   ├── Boundary.cpp
│   ├── Matrix.cpp       # Element stiffness matrices
│   ├── Solver.cpp       # Solve
│   └── Outputs.cpp      # ParaView visualization
├── tests/               #
│   ├── CMakeLists.txt   #
│   └── test.cpp
└── build/               #