#include "../include/Parameter.hpp"
#include <fstream>  // file manipulation
#include <iostream>
#include <format>
#include <utility>

// constructor
ElementInput::ElementInput(std::string name)
  : config_name{std::move(name)}{
    // The variables are already set before this bracket opens.
  }

// display function
void ElementInput::displayElements() const{
    std::cout << std::format("Configuration: {}\n", config_name);
}