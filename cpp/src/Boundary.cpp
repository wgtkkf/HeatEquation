#include "../include/Boundary.hpp"
#include <fstream>  // file manipulation
#include <iostream>
#include <format>
#include <utility>
#include <vector> // for vector
#include <algorithm> // for namespace parseFortranDouble

// constructor
BoundaryInput::BoundaryInput(std::string name)
  : config_name{std::move(name)}{
    // The variables are already set before this bracket opens.
  }

// display function
void BoundaryInput::display_NBC_Numbers() const{
    std::cout << std::format("Configuration: {}\n", config_name);
}

// read rectangle-num.bc
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
    }else if (key == "NumberOfNaturalBoundaryCondition"){
      file >> nbcn;
    }else{
      std::cout << std::format("Warning: unknown parameter '{}' ignored.\n", key);    
    }        
  }  
}

// read rectangle-t.bc
void BoundaryInput::readFixedTemperature(const std::filesystem::path& filepath) {
  std::ifstream file(filepath);

  if (!file.is_open()){
    std::cerr << std::format("Error: file cannot open {}\n", filepath.string());
    return;
  }

  nnbc1.resize(nbc1);  
  temperature.resize(nbc1);  

  // temporal parameters
  int node_id;
  int current_index = 0;
  double t;

  while (file >> node_id >> t){
    nnbc1[current_index] = node_id;    
    temperature[current_index] = t;    
    current_index++;

    if(current_index >= nbc1) break;
  }

  std::cout << std::format("Successfully loaded rectangle-t.bc: {} temperatures.\n", current_index);
}

void BoundaryInput::displayTemperature() const{
  std::cout << std::format("--- Node Temperature ---\n");

  for (size_t i=0; i<nbc1; i++){ // loop with number of nodes
    std::cout << std::format("Node [-] {}: t [K] = {}\n", 
      i + 1, temperature[i]);
  }
}

// read rectangle-n.bc
namespace{
  double parseFortranDouble(std::string token) {
    std::ranges::replace(token, 'd', 'e');
    std::ranges::replace(token, 'D', 'e');
    return std::stod(token);
  }

} // namespace

void BoundaryInput::readNBC(const std::filesystem::path& filepath) {
  std::ifstream file(filepath);

  if (!file.is_open()){
    std::cerr << std::format("Error: file cannot open {}\n", filepath.string());
    return;
  }

  // if you already know the total count
  this->cN.reserve(nbcn);
  this->cM.reserve(nbcn);
  this->cPG.reserve(nbcn);

  // temporal variables
  int n;  
  std::string m, pg;
  int current_index = 0;

  for (int i = 0; i < nbcn; ++i){
    if (file >> n >> m >> pg){
      cN.push_back(n);
      cM.push_back(parseFortranDouble(m));
      cPG.push_back(parseFortranDouble(pg));    
    }
    current_index++;
  }

  std::cout << std::format("Successfully loaded rectangle-n.bc: {}.\n", current_index);
}

void BoundaryInput::displayNBC() const{
  std::cout << std::format("--- N / M / PG ---\n");

  for (size_t i=0; i<nbcn; i++){ // loop with number of nodes
    std::cout << std::format("N [-] {}: M [-] = {}: PG [-] = {}\n", 
      this->cN[i], this->cM[i], this->cPG[i]);
  }
}