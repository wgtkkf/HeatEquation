#include "../include/Parameter.hpp"
#include <iostream>
#include <format>
#include <utility>

// constructor
Input::Input(std::string name)
  : config_name{std::move(name)}{
    // The variables are already set before this bracket opens.
  }

void Input::displayParameters() const{
    std::cout << std::format("Configuration: {}\n", config_name);
}