#pragma once
#include <string>
#include <filesystem> // for filesystem
#include <cstddef> // for std::size_t
#include <vector> // for vector

class BoundaryInput{
  public:       
    std::string config_name{""};  

    // input variables
    int nbc1 = 0; // boundary condition 1
    int nbc2 = 0; // boundary condition 2
    int nbc3 = 0; // boundary condition 3
    int nbcn = 0; // natural boundary condition

    // temperature boundary condition
    std::vector<int> nnbc1;
    std::vector<double> temperature;

    // constructor
    BoundaryInput(std::string name);  

    // read parameters
    void read_NBC_Numbers(const std::filesystem::path& filepath);
    void readFixedTemperature(const std::filesystem::path& filepath);
    void readNBC(const std::filesystem::path& filepath);

    // display
    void display_NBC_Numbers() const; // const guarantees this function only read data.
    void displayTemperature() const;
    
};