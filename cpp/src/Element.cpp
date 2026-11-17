#include "../include/Element.hpp"
#include "../include/Mesh.hpp"
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

// Calculation of the element matrix
void ElementInput::elementMatrix(const MeshInput& mesh) const{
  std::cout << std::format("--- Calculating Element Matrix ---\n");

  for (size_t element=0; element<mesh.ne; ++element){ // loop with number of nodes
  //  std::cout << std::format("Node [-] {}: t [K] = {}\n", 
  //    i + 1, temperature[i]);
  }
}