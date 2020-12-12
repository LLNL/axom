// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <iostream>
#include <string>

#include "axom/inlet.hpp"
#include "axom/slic/core/UnitTestLogger.hpp"

#include "CLI11/CLI11.hpp"

// This example has the following input file:
  /*
  a_simple_bool = true
  a_simple_int = 5,
  a_simple_double = 7.5,
  a_simple_string = 'such simplicity',
  */

bool scalars(const std::string inputFileName)
{
  // Create Inlet Reader that supports Lua input files
  auto lr = std::make_unique<axom::inlet::LuaReader>();

  // Parse example input file
  lr->parseFile(inputFileName);

  // Inlet stores all input file in the Sidre DataStore
  axom::sidre::DataStore ds;

  // Create Inlet with LuaReader and the Sidre Group which Inlet will use
  axom::inlet::Inlet inlet(std::move(lr), ds.getRoot());

  // Define and store the values in the input file
  // _inlet_simple_types_scalar_add_start
  inlet.addBool("a_simple_bool", "A description of a_simple_bool");
  inlet.addInt("a_simple_int", "A description of a_simple_int");
  inlet.addDouble("a_simple_double", "A description of a_simple_double").required();
  inlet.addString("a_simple_string", "A description of a_simple_string");

  inlet.addInt("a_defaulted_int", "An int that has a default value").defaultValue(17);

  inlet.addString("does_not_exist", "Should not be in your input file");
  // _inlet_simple_types_scalar_add_end


  // Access values stored in the Datastore via Inlet
  // _inlet_simple_types_scalar_access_start
  bool a_simple_bool = inlet["a_simple_bool"];
  int a_simple_int = inlet["a_simple_int"];
  double a_simple_double = inlet["a_simple_double"];
  std::string a_simple_string = inlet["a_simple_string"];

  int a_defaulted_int = inlet["a_defaulted_int"];
  // _inlet_simple_types_scalar_access_end

  // _inlet_simple_types_scalar_contains_end
  std::string does_not_exist = "Yay! Was not present!";
  if (inlet.contains("does_not_exist"))
  {
    does_not_exist = "Uh oh! This should not happen!";
  }
  // _inlet_simple_types_scalar_contains_end

  std::cout << "a_simple_bool = " << a_simple_bool << std::endl;
  std::cout << "a_simple_int = " << a_simple_int << std::endl;
  std::cout << "a_simple_double = " << a_simple_double << std::endl;
  std::cout << "a_simple_string = " << a_simple_string << std::endl;
  std::cout << "a_defaulted_int = " << a_defaulted_int << std::endl;
  std::cout << "does_not_exist = " << does_not_exist << std::endl;

  return true;
}


int main()
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  // This is a generic basic SLIC logger
  axom::slic::UnitTestLogger logger;


  // Handle command line arguments
  CLI::App app {"Basic example of Axom's Inlet component with simple types"};
  std::string inputFileName;
  auto opt = app.add_option("--file", inputFileName, "Path to input file");
  opt->check(CLI::ExistingFile);
  CLI11_PARSE(app, argc, argv);

  // Small example of basic scalar types
  if(!scalars(inputFileName))
  {
    return 1;
  }

  return 0;
}
