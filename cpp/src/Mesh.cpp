#include "../include/Mesh.hpp"
#include <fstream>  // file manipulation
#include <iostream>
#include <format>
#include <utility>
#include <string> // for char variable

// constructor
MeshInput::MeshInput(std::string name)
  : config_name{std::move(name)}{
    // The variables are already set before this bracket opens.
  }

void MeshInput::displayParameters() const{
    std::cout << std::format("Configuration: {}\n", config_name);
}

// read rectangle1.msh
void MeshInput::readMeshParameters(const std::filesystem::path& filepath) {
  std::ifstream file(filepath);

  if (!file.is_open()){
    std::cerr << std::format("Error: file cannot open {}\n", filepath.string());
    return;
  }

  std::string key;
  // read file by keyword
  while (file >> key){
    if(key=="GeometryShape"){
      file >> geometry;
    }else if (key == "Nodes"){
      file >> np;
    }else if (key == "Elements"){
      file >> ne;
    }else if (key == "TimeStep"){
      file >> dt;
    }else{
      std::cout << std::format("Warning: unknown parameter '{}' ignored.\n", key);    
    }
  } 

  std::cout << std::format("Successfully loaded rectanble1.msh: node, element, time step.\n");
}

// read rectangle2.msh
void MeshInput::readMeshCoordinates(const std::filesystem::path& filepath) {
  std::ifstream file(filepath);

  if (!file.is_open()){
    std::cerr << std::format("Error: file cannot open {}\n", filepath.string());
    return;
  }

  x_coords.resize(np);
  y_coords.resize(np);

  // temporal parameters
  int node_id;
  int current_index = 0;
  double x, y;

  while (file >> node_id >> x >> y){
    x_coords[current_index] = x;
    y_coords[current_index] = y;
    current_index++;

    if(current_index >= np) break;
  }

  std::cout << std::format("Successfully loaded rectanble2.msh: x & y coordinates.\n", current_index);
}

void MeshInput::displayCoordinates(int total_nodes) const{
  std::cout << std::format("--- Node coordinates ---\n");

  for (size_t i=0; i<total_nodes; i++){ // loop with number of nodes
    std::cout << std::format("Node {}: x = {}, y = {}\n", 
      i + 1, x_coords[i], y_coords[i]);
  }
}

// read rectangle3.msh
void MeshInput::readMeshNodes(const std::filesystem::path& filepath) {
  std::ifstream file(filepath);

  if (!file.is_open()){
    std::cerr << std::format("Error: file cannot open {}\n", filepath.string());
    return;
  }

  nop.resize(ne);

  // temporal parameters
  int element_id;
  int current_index = 0;
  double a, b, c;

  while (file >> element_id >> a >> b >> c){
    nop[current_index][0] = a;
    nop[current_index][1] = b;
    nop[current_index][2] = c;
        
    current_index++;

    if(current_index >= ne) break;
  }

<<<<<<< HEAD
  std::cout << std::format("Successfully loaded rectanble3.msh: node information.\n", current_index);
}

void MeshInput::displayNodes(int total_elements) const{
  std::cout << std::format("--- Nodes of each element ---\n");

  for (size_t i=0; i<total_elements; i++){ // loop with number of element
    std::cout << std::format("Element {}: a = {}, b = {}, c = {}\n", 
      i + 1, nop[i][0], nop[i][1], nop[i][2]);
  }
=======
  std::cout << std::format("Successfully loaded x & ycoordinates.\n", current_index);
>>>>>>> f3e1f5e2ab401f2575ee0661bd2ad6606cfb7075
}