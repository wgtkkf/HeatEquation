#include <gtest/gtest.h>
#include <fstream> // file manipulation
#include <cstdio>
#include <string>
#include "../include/Mesh.hpp"

// 1. Define a Fixture class
class MeshInputTest : public ::testing::Test{
  protected:
    std::string test_file = "mock_rectangle1.msh";

    void SetUp() override {
      std::ofstream out(test_file);
      out << "GeometryShape Square\n";
      out << "Nodes 1\n";
      out << "Elements 2\n";
      out << "TimeStep 0.1\n";
      out << "UnDefinedKeyword 999\n";
      out.close();
    }

    // TearDown runs automatically aafter every TEST_F, even if it crashes
    void TearDown() override {
      std::remove(test_file.c_str());
    }
  };


TEST_F(MeshInputTest, ReadValidMeshParameters) {
  
  // 1. Execution: run moduls to be tested
  MeshInput input("MockParameters");
  input.readMeshParameters(test_file);

  // 2. Verification: assert the variables were populated correctly
  EXPECT_EQ(input.geometry, "Square");
  EXPECT_EQ(input.getNumNodes(), 1);
  EXPECT_EQ(input.getNumElements(), 2);
  EXPECT_DOUBLE_EQ(input.getDt(), 0.1);    
}

// Test Case 2: missing file, error handling
TEST_F(MeshInputTest, ReadMissingFile){
  MeshInput input("MockDummyParameters");

  // read a file that does not exist
  input.readMeshParameters("non_existent_file.msh")  ;

  EXPECT_EQ(input.getNumNodes(), 0);
}