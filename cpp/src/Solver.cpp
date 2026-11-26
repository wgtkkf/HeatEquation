#include "../include/Solver.hpp"
#include "../include/Mesh.hpp"

#include <cstddef>
#include <format>
#include <iostream>
#include <stdexcept>
#include <utility>

Solver::Solver(std::string name)
    : config_name{std::move(name)} {}

void Solver::displaySolver() const {
    std::cout << std::format("Configuration: {}\n", config_name);
}

void Solver::initialize(int num_nodes, double initial_temp) {
    if (num_nodes <= 0) {
        throw std::invalid_argument(
            std::format("num_nodes must be positive, got {}", num_nodes));
    }
    T_ = Eigen::VectorXd::Constant(num_nodes, initial_temp);
}

void Solver::advanceStep([[maybe_unused]] const MeshInput& mesh) {
    // TODO: Assemble global matrices, apply boundary conditions, and solve for T_
}

std::span<const double> Solver::getTemperature() const noexcept {
    return {T_.data(), static_cast<std::size_t>(T_.size())};
}