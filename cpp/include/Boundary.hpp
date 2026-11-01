#pragma once
#include <string>
#include <filesystem> // for filesystem
#include <cstddef> // for std::size_t

class BoundaryInput{
  public:       
    std::string config_name{""};

    // input variables
    int nbc1 = 0; // natural boundary condition 1
    int nbc2 = 0; // natural boundary condition 2
    int nbc3 = 0; // natural boundary condition 3
    

    // constructor
    BoundaryInput(std::string name);  

    // display
    void display_NBC_Numbers() const; // const guarantees this function only read data.

    // read parameters
    void read_NBC_Numbers(const std::filesystem::path& filepath);
};