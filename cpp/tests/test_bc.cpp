#include <gtest/gtest.h>
#include <fstream> // file manipulation
#include <cstdio>
#include <string>
#include "../include/Boundary.hpp"

// 1. Define a Fixture class
class MeshInputTest : public ::testing::Test{
  protected:
    std::string test_file = "mock_rectangle-num.bc";

    void SetUp() override {
      std::ofstream out(test_file);
      out << "GeometryShape Square\n";
      out << "NaturalBoundaryConditionNumber1 100\n";
      out << "NaturalBoundaryConditionNumber2 50\n";
      out << "NaturalBoundaryConditionNumber3 10\n";
      out.close();
    }

    // TearDown runs automatically aafter every TEST_F, even if it crashes
    void TearDown() override {
      std::remove(test_file.c_str());
    }
  };

TEST_F(MeshInputTest, ReadValidMeshParameters) {
  
  // 1. Execution: run moduls to be tested
  BoundaryInput input("MockParameters");
  input.read_NBC_Numbers(test_file);

  // 2. Verification: assert the variables were populated correctly
  EXPECT_EQ(input.geometry, "Square");
  EXPECT_EQ(input.nbc1, 100);
  EXPECT_EQ(input.nbc2, 50);
  EXPECT_EQ(input.nbc3, 10);
}

// Test Case 2: missing file, error handling
TEST_F(BoundaryinputTest, ReadMissingFile){
  BoundaryInput input("MockDummyParameters");

  // Set a default dummy values to prove they don't get overwritten by garbage data
  input.nbc1 = 0;

  // read a file that does not exist
  input.read_NBC_Numbers("non_existent_file.msh")  ;

  EXPECT_EQ(input.nbc1, 0);
}