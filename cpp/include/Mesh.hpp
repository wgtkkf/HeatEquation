#pragma once
#include <string>
#include <filesystem> // for filesystem
#include <cstddef> // for std::size_t
#include <string> // for char variable
#include <vector> // for vector

class MeshInput{
  private:
    static constexpr std::size_t size0 = 3;
    static constexpr std::size_t size1 = 2;    

    // input variables
    std::string geometry{};

  public:    
    // input variables
    double np = 0.0; /* nodes */
    double ne = 0.0; /* elements */
    double dt = 0.0; /* time step */

    // mesh coordination
    std::vector<double> x_coords;
    std::vector<double> y_coords;

    std::string config_name{""};

    // constructor
    MeshInput(std::string name);
    // display
    void displayParameters() const; // const guarantees this function only read data.

    // read parameters
    void readMeshParameters(const std::filesystem::path& filepath);

    // read x & y coordinates
    void readMeshCoordinates(const std::filesystem::path& filepath);

    void displayCoordinates(int total_nodes) const; // The 'const' guarantees this method not altering data
};