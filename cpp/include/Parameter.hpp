#pragma once
#include <string>

class Input{
  public:
    int size0 = 3;
    std::string config_name{""};

    // constructor
    Input(std::string name);

    // display
    void displayParameters() const; // const guarantees this function only read data.
};