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

// read parameter.txt
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
}