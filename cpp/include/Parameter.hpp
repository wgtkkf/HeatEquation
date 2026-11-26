#pragma once
#include <string>
#include <filesystem> // for filesystem
#include <cstddef> // for std::size_t

class ParameterInput{
  public:       
    std::string config_name{""};

    // input variables
    double thermal_conductivity = 0.0;
    double density = 0.0;
    double specific_heat = 0.0;

    [[nodiscard]] double getC1() const noexcept { return c1; }

    // constructor
    ParameterInput(std::string name);  

    // display
    void displayParameters() const; // const guarantees this function only read data.

    // read parameters
    void readParameters(const std::filesystem::path& filepath);

  private:
    double c1{0.0};
};