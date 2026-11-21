#pragma once
#include <string>
#include <filesystem> // for filesystem
#include <cstddef> // for std::size_t
#include "Mesh.hpp" // for elementMatrix function

class ElementInput{
  public:       
    std::string config_name{""};    
    ElementInput(std::string name); // constructor    

    //        
    using ElementMat = Eigen::Matrix<double, MeshInput::size0, MeshInput::size0>;
    std::vector<Eigen::MatrixXd> ass; // One 2D matrix per element    
    Eigen::MatrixXd emm;    
    Eigen::VectorXd lmm;
    Eigen::VectorXd ilmm;

    std::vector<std::array<double, MeshInput::size0>> eb;
    std::vector<std::array<double, MeshInput::size0>> ec;

    //
    void displayElements() const;
    void elementMatrix(const MeshInput& mesh);
};