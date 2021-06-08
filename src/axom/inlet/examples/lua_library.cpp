// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>

#include "axom/inlet.hpp"
#include "axom/slic/core/SimpleLogger.hpp"

int main()
{
  const std::string test_file_name = "load_library_test_file";

  const std::string input = R"(
  read_str = function ()
      file = io.open('load_library_test_file', 'r')
      return file:read()
  end
  )";

  // Inlet requires a SLIC logger to be initialized to output runtime information
  // This is a generic basic SLIC logger
  axom::slic::SimpleLogger logger;

  // Write test file
  std::ofstream myfile;
  myfile.open(test_file_name);
  myfile << "test_string";
  myfile.close();

  // _inlet_io_library_add_start
  // Create Inlet Reader that supports Lua input files
  auto lr = std::make_unique<axom::inlet::LuaReader>();

  // Load extra io Lua library
  lr->solState().open_libraries(sol::lib::io);

  // Parse example input string
  lr->parseString(input);
  // _inlet_io_library_add_end

  // Inlet stores all input information in the Sidre DataStore
  axom::sidre::DataStore ds;

  // Create Inlet with LuaReader and the Sidre Group which Inlet will use
  axom::inlet::Inlet myinlet(std::move(lr), ds.getRoot());

  // Define and store the values in the input file
  myinlet.addFunction("read_str",
                      axom::inlet::FunctionTag::String,  // Return type
                      {},                                // Argument types
                      "The function reads a double from a file");

  auto result = myinlet["read_str"].call<std::string>();

  if(result != "test_string")
  {
    std::cerr << "Failed to read 'test_string' from test file." << std::endl;
    return 1;
  }
  std::cout << "Successfully read 'test_string' from test file." << std::endl;

  // Clean up test file
  remove(test_file_name.c_str());

  return 0;
}
