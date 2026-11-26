#include "../include/Output.hpp"
#include <fstream>  // file manipulation
#include <iostream>
#include <format>
#include <utility>

// constructor
OutputResult::OutputResult(std::string name)
  : config_name{std::move(name)}{
    // The variables are already set before this bracket opens.
  }

// display function
void OutputResult::OutputResult() const{
    std::cout << std::format("Configuration: {}\n", config_name);
}