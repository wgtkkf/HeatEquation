#pragma once
#include <string>

#include <cstddef> // for std::size_t

class MeshInput{
  private:
    static constexpr std::size_t size0 = 3;
    static constexpr std::size_t size1 = 2;    

  public:
    std::string config_name{""};

    // constructor
    MeshInput(std::string name);
    // display
    void displayParameters() const; // const guarantees this function only read data.
};