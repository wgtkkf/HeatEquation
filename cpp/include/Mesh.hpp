#pragma once
#include <string>
#include <filesystem> // for filesystem
#include <cstddef> // for std::size_t
#include <string> // for char variable
#include <vector> // for vector

class MeshInput{
  private:    
    // input variables
    std::string geometry{};

  public: 
    static constexpr std::size_t size0 = 3;
    static constexpr std::size_t size1 = 2;       

    // input variables
    double np = 0.0; /* nodes */
    double ne = 0.0; /* elements */
    double dt = 0.0; /* time step */

    // mesh coordination
    std::vector<double> x_coords;
    std::vector<double> y_coords;

    // nodes on each element
    std::vector<double> nop(ne, size0, 0.0);

    std::string config_name{""};

    // constructor
    MeshInput(std::string name);

    // display
    void displayParameters() const; // const guarantees this function only read data.

    // read parameters
    void readMeshParameters(const std::filesystem::path& filepath);

    // read x & y coordinates
    void readMeshCoordinates(const std::filesystem::path& filepath);
    
    // show x & y coordinates
    void displayCoordinates(int total_nodes) const; // The 'const' guarantees this method not altering data

    // read node information
    void readMeshNodes(const std::filesystem::path& filepath);
};