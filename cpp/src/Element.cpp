#include "../include/Element.hpp"
#include "../include/Mesh.hpp"
#include <fstream>  // file manipulation
#include <iostream>
#include <format>
#include <utility>

#include <format>
#include <stdexcept> // for runtime_error

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
void ElementInput::elementMatrix(const MeshInput& mesh) {
  std::cout << std::format("--- Calculating Element Matrix ---\n");

  /* initialization */
  emm.setZero(mesh.ne, mesh.size0);
  ass.assign(mesh.ne, Eigen::MatrixXd::Zero(mesh.size0, mesh.size0));  

  for (size_t element=0; element<mesh.ne; ++element){ // loop with number of nodes
    /* calculation area of element */
    /* n1 -- number of 1st node */
		/* n2 -- number of 2nd node */
		/* n3 -- number of 3rd node */
    
    const auto [n1, n2, n3] = mesh.nop[element]; // Unpack all 3 nodes in one line (replaces n1 = nop[elem][0], etc.)
    // The above is equivalent with the below three lines
    // const int n1 = mesh.nop[elem][0];
    // const int n2 = mesh.nop[elem][1];
    // const int n3 = mesh.nop[elem][2];

    const double x1 = mesh.x_coords.at(n1);
    const double y1 = mesh.y_coords.at(n1);
    const double x2 = mesh.x_coords.at(n2);
    const double y2 = mesh.y_coords.at(n2);
    const double x3 = mesh.x_coords.at(n3);
    const double y3 = mesh.y_coords.at(n3);

    // Inside your element loop in Element.cpp:
    const double area2 = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
    const double area = 0.5 * area2;    
    std::cout << std::format("elem = {}, area = {:.6f}\n", element, area);

    if (area <= 0.0) {
      throw std::runtime_error(
        std::format("Error: element {} has non-positive area ({:.6e})", element, area));
      }
  }
}