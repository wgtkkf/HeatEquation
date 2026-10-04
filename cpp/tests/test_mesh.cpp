#include <gtest/gtest.h>
#include <fstream> // file manipulation
#include <cstdio>
#include <string>
#include "../include/Mesh.hpp"

TEST(MeshInputTest, ReadValidMeshParameters) {
  // 1. Setup: create a temporary mock file with expected keys
  std::string test_file = "mock_rectangle1.msh";
  std::ofstream out(test_file);
  out << "GeometryShape Square\n";
  out << "Nodes 1\n";
  out << "Elements 2\n";
  out << "TimeStep 0.1\n";
  out << "UnDefinedKeyword 999\n";
  out.close();

  // 2. Execution: run moduls to be tested
  MeshInput input("MockParameters");
  input.readMeshParameters(test_file);

  // 3. Verification: assert the variables were populated correctly
  EXPECT_EQ(input.geometry, "Square");
  EXPECT_EQ(input.np, 1);
  EXPECT_EQ(input.ne, 2);
  EXPECT_DOUBLE_EQ(input.dt, 0.1);
  
  // 4. clean up the file
  std::remove(test_file.c_str());
}

// Test Case 2: missing file, error handling
TEST(MeshInputTest, ReadMissingFile){
  MeshInput input("MockDummyParameters");

  // Set a default dummy values to prove they don't get overwritten by garbage data
  input.np = 0;

  // read a file that does not exist
  input.readMeshParameters("non_existent_file.msh")  ;

  EXPECT_EQ(input.np, 0);
}