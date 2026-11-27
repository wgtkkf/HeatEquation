#include "../include/Output.hpp"
#include "../include/Mesh.hpp"
#include <fstream>  // file manipulation
#include <iostream>
#include <format>
#include <utility>
#include <stdexcept>

// constructor
OutputResult::OutputResult(std::string name,
                           std::filesystem::path output_dir,
                           std::string file_prefix)
  : config_name{std::move(name)},
  dir_{std::move(output_dir)},
  prefix_{std::move(file_prefix)}{}

// display function
void OutputResult::displayOutput() const{
    std::cout << std::format("Configuration: {}\n", config_name)
              << std::format("Output directory: {}\n", dir_.string())
              << std::format("File prefix: {}\n", prefix_);
}

void OutputResult::writeVtk(int stepnum,
                            const MeshInput& mesh,
                            std::span<const double> temperature) const {
    /* initialization of a parameter */
    constexpr int kVtkTriangle{5};

    /* consistency of np, x_coords, and y_coords */
    const auto& x = mesh.getXCoords();
    const auto& y = mesh.getYCoords();
    const auto np = mesh.getNumNodes();

    /* new description for C++ version, not C version */
    const auto& conn = mesh.getConnectivity();
    const auto ne = mesh.getNumElements();

    /* size consistency check */
    if (std::cmp_not_equal(np, x.size()) || std::cmp_not_equal(np, y.size())) {
      throw std::runtime_error(std::format(
        "Mesh node count mismatch: np={}, x_coords.size()={}, y_coords.size()={}",
        np, x.size(), y.size()));
    }

    if (std::cmp_not_equal(ne, conn.size())) {
        throw std::runtime_error(std::format(
            "Mesh element count mismatch: ne={}, connectivity.size()={}",
            ne, conn.size()));
    }

    if (std::cmp_not_equal(np, temperature.size())) {
        throw std::runtime_error(std::format(
            "Temperature field size mismatch: np={}, temperature.size()={}",
            np, temperature.size()));
    }

    if (!dir_.empty()) {
        std::filesystem::create_directories(dir_);
    }

    /* 7.1 creation of a file */
    const std::filesystem::path filename =
        dir_ / std::format("{}{}.vtk", prefix_, stepnum);
    
    std::ofstream fp(filename);
    if (!fp) {
      throw std::runtime_error(
        std::format("Error: cannot open file {}", filename.string()));
      }
    
    /* 7.2 write wods as VTK header of PARAVIEW format */
    fp << "# vtk DataFile Version 2.0\n"
       << "temperature\n"
       << "ASCII\n"
       << "DATASET UNSTRUCTURED_GRID\n"
       << std::format("POINTS {} double\n", np);

    // coord information
    for (std::size_t i = 0; i < static_cast<std::size_t>(np); ++i) {
      fp << std::format("{:.6f} {:.6f} 0.0\n", x[i], y[i]);
    }

    // element informat
    const auto total_cell_entries = static_cast<std::size_t>(ne) * (MeshInput::size0 + 1);
    fp << std::format("CELLS {} {}\n", ne, total_cell_entries);

    // this line may need double check
    for (const auto& elem : conn) {
      fp << std::format("{} {} {} {}\n", MeshInput::size0, elem[0], elem[1], elem[2]);
    }

    fp << std::format("CELL_TYPES {}\n", ne);
    for (int i = 0; i < ne; ++i) {
      fp << std::format("{}\n", kVtkTriangle);
    }

    /* temperature data */
    fp << std::format("POINT_DATA {}\n", np)
       << "SCALARS temperature double 1\n"
       << "LOOKUP_TABLE default\n";

    for (const double temp : temperature) {
        fp << std::format("{:.6f}\n", temp);
    }

    fp << std::format("CELL_DATA {}\n", ne);    
}