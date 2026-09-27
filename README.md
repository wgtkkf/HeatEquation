## About this repository
Two-dimensional finite element modeling of heat equasion.

<<<<<<< HEAD
=======
<<<<<<< Updated upstream
# Pre-requisits
1. Create two folders: vtk and excel
2. Run: gcc source.c -o source -lm
=======
>>>>>>> dev_1
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
<<<<<<< HEAD
│   ├── source.cpp         # entry point
=======
│   ├── source.cpp       # entry point
>>>>>>> dev_1
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
<<<<<<< HEAD
=======
├── inputs/              # external input files
│   ├── parameters.txt   # 
>>>>>>> dev_1
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

Remove all the files in your build folder if the cmake command does not work.
```text
cd build
rm -rf *
```

### CMake & Docker
```text
docker build -t simulation .
docker run --rm simulation
<<<<<<< HEAD
```
=======
```
>>>>>>> Stashed changes
>>>>>>> dev_1
