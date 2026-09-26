#pragma once
#include <string>

#include <cstddef> // for std::size_t

class Input{
  public:
    static constexpr std::size_t size0 = 3;
    static constexpr std::size_t size1 = 2;
    
    std::string config_name{""};

    // constructor
    Input(std::string name);

    // display
    void displayParameters() const; // const guarantees this function only read data.
};