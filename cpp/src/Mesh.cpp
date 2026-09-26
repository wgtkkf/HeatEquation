#include "../include/Mesh.hpp"
#include <iostream>
#include <format>
#include <utility>

// constructor
MeshInput::MeshInput(std::string name)
  : config_name{std::move(name)}{
    // The variables are already set before this bracket opens.
  }

void MeshInput::displayParameters() const{
    std::cout << std::format("Configuration: {}\n", config_name);
}