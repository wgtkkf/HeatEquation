## About this repository
Two-dimensional finite element modeling of heat equasion.

### A prerequisit for C-language code
Please create two folders: vtk and excel

### Run C-language code
1 gcc source.c -o source -lm  
2 ./source

## Boundary conditions
Fixed temperature, a schematic is in preparation.

### Analytical solution
An analytical solution was given by a python script.

## C++ version
```text
cpp_version/
├── CMakeLists.txt       # 
├── Dockerfile           # 
├── include/             # .hpp files
│   ├── Timer.hpp
│   ├── Parameter.hpp
│   ├── Mesh.hpp
│   ├── Initial.hpp
│   ├── Boundary.hpp
│   ├── Matrix.hpp
│   ├── Solver.hpp      
│   └── Outputs.hpp      
├── src/                 # .cpp files
│   ├── main.cpp         # entry point
│   ├── Timer.cpp
│   ├── Parameter.cpp
│   ├── Mesh.cpp
│   ├── Initial.cpp
│   ├── Boundary.cpp
│   ├── Matrix.cpp       # Element stiffness matrices
│   ├── Solver.cpp       # Solve
│   └── Outputs.cpp      # ParaView visualization
├── tests/               #
│   ├── test1.txt        #
│   └── test2.cpp
└── build/               #
```

### CMake
```text
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=g++-13 ..
make
./bmt
```

If the cmake command does not work, remove all the files in your build folder.
```text
rm -rf *
```

### CMake & Docker
```text
docker build -t simulation .
docker run --rm simulation
```