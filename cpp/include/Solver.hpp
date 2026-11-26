#pragma once

#include <Eigen/Dense>
#include <span>
#include <string>

class MeshInput;

class Solver {
public:
    explicit Solver(std::string name);

    void displaySolver() const;
    void initialize(int num_nodes, double initial_temp = 0.0);
    void advanceStep(const MeshInput& mesh);

    [[nodiscard]] std::span<const double> getTemperature() const noexcept;

private:
    std::string config_name;
    Eigen::VectorXd T_;
};