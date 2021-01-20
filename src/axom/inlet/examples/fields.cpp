// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <iostream>
#include <string>

#include "axom/inlet.hpp"
#include "axom/slic/core/SimpleLogger.hpp"

#include "CLI11/CLI11.hpp"

/* Input file snippet used for documentation
//_inlet_simple_types_fields_input_start

a_simple_bool = true
a_simple_int = 5
a_simple_double = 7.5
a_simple_string = 'such simplicity'

//_inlet_simple_types_fields_input_end
*/

int main()
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  // This is a generic basic SLIC logger
  axom::slic::SimpleLogger logger;

  std::string input = R""""(
    a_simple_bool = true
    a_simple_int = 5
    a_simple_double = 7.5
    a_simple_string = 'such simplicity'
  )"""";

  // Create Inlet Reader that supports Lua input files
  auto lr = std::make_unique<axom::inlet::LuaReader>();

  // Parse example input file string
  // Note: the Reader class also supports parseFile(std::string filepath)
  lr->parseString(input);

  // Inlet stores all input file in the Sidre DataStore
  axom::sidre::DataStore ds;

  // Create Inlet with LuaReader and the Sidre Group which Inlet will use
  axom::inlet::Inlet inlet(std::move(lr), ds.getRoot());

  // Define and store the values in the input file
  // _inlet_simple_types_fields_add_start
  inlet.addBool("a_simple_bool", "A description of a_simple_bool");
  inlet.addInt("a_simple_int", "A description of a_simple_int");
  inlet.addDouble("a_simple_double", "A description of a_simple_double").required();
  inlet.addString("a_simple_string", "A description of a_simple_string");

  inlet.addInt("a_defaulted_int", "An int that has a default value").defaultValue(17);

  inlet.addString("does_not_exist", "Should not be in your input file");
  // _inlet_simple_types_fields_add_end

  // Access values stored in the Datastore via Inlet
  // _inlet_simple_types_fields_access_start
  bool a_simple_bool = inlet["a_simple_bool"];
  int a_simple_int = inlet["a_simple_int"];
  double a_simple_double = inlet["a_simple_double"];
  std::string a_simple_string = inlet["a_simple_string"];

  int a_defaulted_int = inlet["a_defaulted_int"];
  // _inlet_simple_types_fields_access_end

  std::cout << "a_simple_bool = " << a_simple_bool << std::endl;
  std::cout << "a_simple_int = " << a_simple_int << std::endl;
  std::cout << "a_simple_double = " << a_simple_double << std::endl;
  std::cout << "a_simple_string = " << a_simple_string << std::endl;
  std::cout << "a_defaulted_int = " << a_defaulted_int << std::endl;

  // _inlet_simple_types_fields_contains_start
  if(inlet.contains("does_not_exist"))
  {
    std::cout << "Error: Inlet should not have contained key 'does_not_exist' "
              << std::endl;
    return 1;
  }
  else
  {
    std::cout << "Key 'does_not_exist' did not exist in Inlet. Success!"
              << std::endl;
  }
  // _inlet_simple_types_fields_contains_end

  return 0;
}
