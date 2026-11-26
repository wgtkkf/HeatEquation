#pragma once
#include <string>
#include <filesystem> // for filesystem
#include <cstddef> // for std::size_t
#include <string> // for char variable
#include <vector> // for vector
#include <Eigen/Dense> // for matrix operation

class MeshInput{

  public:
    // input variables
    std::string geometry{};

    static constexpr std::size_t size0{3};
    static constexpr std::size_t size1{2};    

    [[nodiscard]] int getNumNodes() const noexcept { return np; }    /* nodes */
    [[nodiscard]] int getNumElements() const noexcept { return ne; } /* elements */
    [[nodiscard]] double getDt() const noexcept { return dt; }       /* time step */

    // mesh coordination    
    [[nodiscard]] const std::vector<double>& getXCoords() const noexcept { return x_coords; }
    [[nodiscard]] const std::vector<double>& getYCoords() const noexcept { return y_coords; }
    
    //  node information of each element
    std::vector<std::array<int, size0>> nop;

    std::string config_name{""};

    // Modules
    // constructor
    MeshInput(std::string name);

    // display
    void displayParameters() const; // const guarantees this function only read data.

    // read parameters
    void readMeshParameters(const std::filesystem::path& filepath);

    // read x & y coordinates
    void readMeshCoordinates(const std::filesystem::path& filepath);
    void displayCoordinates(int total_nodes) const; // The 'const' guarantees this method not altering data
    
    // read nodes on each elements
    void readMeshNodes(const std::filesystem::path& filepath);
    void displayNodes(int total_nodes) const; // The 'const' guarantees this method not altering data

    // This is used in Output.cpp
    [[nodiscard]] const auto& getConnectivity() const noexcept {
    return nop; }

  private:
    // input variables
    int np{0};      /* nodes */
    int ne{0};      /* elements */
    double dt{0.0}; /* time step */

    std::vector<double> x_coords;
    std::vector<double> y_coords;
};