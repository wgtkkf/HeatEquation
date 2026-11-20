#pragma once
#include <string>
#include <filesystem> // for filesystem
#include <cstddef> // for std::size_t
#include "Mesh.hpp" // for elementMatrix function

class ElementInput{
  public:       
    std::string config_name{""};    
    ElementInput(std::string name); // constructor

    void displayElements() const;
    void elementMatrix(const MeshInput& mesh);

    //
    Eigen::MatrixXd emm;
    std::vector<Eigen::MatrixXd> ass; // One 2D matrix per element    
};