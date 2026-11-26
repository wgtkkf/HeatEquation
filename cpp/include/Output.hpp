#pragma once
#include <string>
#include <filesystem> // for filesystem
#include <cstddef>    // for std::size_t
#include <utility>
#include <format>

class MeshInput;

class OutputResult{
  public:           
    OutputResult(
      std::string name,
      std::filesystem::path output_dir,
      std::string file_prefix);        

    void display() const;
    void writeVtk(int stepnum, const MeshInput& mesh) const;

private:
    std::string config_name;
    std::filesystem::path dir_;
    std::string prefix_;

};