#include "../include/Boundary.hpp"
#include <fstream>  // file manipulation
#include <iostream>
#include <format>
#include <utility>

// constructor
BoundaryInput::BoundaryInput(std::string name)
  : config_name{std::move(name)}{
    // The variables are already set before this bracket opens.
  }

// display function
void BoundaryInput::display_NBC_Numbers() const{
    std::cout << std::format("Configuration: {}\n", config_name);
}

// read parameter.txt
void BoundaryInput::read_NBC_Numbers(const std::filesystem::path& filepath) {
  std::ifstream file(filepath);

  if (!file.is_open()){
    std::cerr << std::format("Error: file cannot open {}\n", filepath.string());
    return;
  }

  std::string key;

  // read file by keyword
  while (file >> key){
    if(key=="NumberOfBoundaryCondition1"){
      file >> nbc1;
    }else if (key == "NumberOfBoundaryCondition2"){
      file >> nbc2;
    }else if (key == "NumberOfBoundaryCondition3"){
      file >> nbc3;
    }else{
      std::cout << std::format("Warning: unknown parameter '{}' ignored.\n", key);    
    }        
  }  
}