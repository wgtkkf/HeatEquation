#pragma once
#include <filesystem> // for filesystem
#include <span>
#include <string>

class MeshInput;

class OutputResult{
  public:           
    OutputResult(
      std::string name,
      std::filesystem::path output_dir,
      std::string file_prefix);        

    void displayOutput() const;
    void writeVtk(int stepnum,
                  const MeshInput& mesh,
                  std::span<const double> temperature) const;

private:
    std::string config_name;
    std::filesystem::path dir_;
    std::string prefix_;
};