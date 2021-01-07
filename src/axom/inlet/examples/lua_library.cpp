// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>

#include "axom/inlet.hpp"
#include "axom/slic/core/UnitTestLogger.hpp"

#include "CLI11/CLI11.hpp"

const char* TEST_FILE_NAME = "load_library_test_file";

// TODO: add string and None to functions?
const char* LUA_STRING = R"PREFIX(
read_str = function (x,y,z)
    file = io.open('load_library_test_file', 'r')
    return tonumber(file:read())
end
)PREFIX";

bool load_library()
{
  // _inlet_io_library_add_start
  // Create Inlet Reader that supports Lua input files
  auto lr = std::make_unique<axom::inlet::LuaReader>();

  // Load extra io Lua library
  lr->solState().open_libraries(sol::lib::io);

  // Parse example input file
  lr->parseString(LUA_STRING);
  // _inlet_io_library_add_end

  // Inlet stores all input file in the Sidre DataStore
  axom::sidre::DataStore ds;

  // Create Inlet with LuaReader and the Sidre Group which Inlet will use
  axom::inlet::Inlet myinlet(std::move(lr), ds.getRoot());

  // Define and store the values in the input file
  myinlet.addFunction("read_str",
                      axom::inlet::FunctionType::Double,   // Return type
                      {axom::inlet::FunctionType::Vec3D}, // Argument types (unused)
                      "The function reads a double from a file");

  // call the function with an unimportant number
  // TODO: support no arguments?
  //double result = myinlet["coef"].call<double>(5.0);
  double result = myinlet["read_str"].call<double>(axom::primal::Vector3D{3, 5, 7});

  if (result != 7.4) {
    std::cerr << "Failed to read '7.4' from test file." << std::endl;  
  }
  std::cout << "Successfully read '7.4' from test file." << std::endl;
  return true;
}


int main()
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  // This is a generic basic SLIC logger
  axom::slic::UnitTestLogger logger;

  // Write test file
  std::ofstream myfile;
  myfile.open(TEST_FILE_NAME);
  myfile << "7.4";
  myfile.close();

  // Small example of loading an extra Lua library
  load_library();

  // Clean up test file
  if (!remove(TEST_FILE_NAME)) {
    return 1;
  }

  return 0;
}
