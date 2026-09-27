#include "../include/Parameter.hpp"
#include <fstream>  // file manipulation
#include <iostream>
#include <format>
#include <utility>

// constructor
ParameterInput::ParameterInput(std::string name)
  : config_name{std::move(name)}{
    // The variables are already set before this bracket opens.
  }

// display function
void ParameterInput::displayParameters() const{
    std::cout << std::format("Configuration: {}\n", config_name);
}

// read parameter.txt
void ParameterInput::readParameters(const std::filesystem::path& filepath) {
  std::ifstream file(filepath);

  if (!file.is_open()){
    std::cerr << std::format("Error: file cannot open {}\n", filepath.string());
    return;
  }

  std::string key;

  // read file by keyword
  while (file >> key){
    if(key=="ThermalConductivity"){
      file >> thermal_conductivity;
    }else if (key == "Density"){
      file >> density;
    }else if (key == "SpecificHeat"){
      file >> specific_heat;
    }else{
      std::cout << std::format("Warning: unknown parameter '{}' ignored.\n", key);    
    }
  }
}