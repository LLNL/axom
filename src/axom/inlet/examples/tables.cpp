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
//_inlet_simple_types_tables_input_start

driver = {
    name = "Speed Racer",
    car = {
        make = "BestCompany",
        seats = 2,
        horsepower = 200
    }
}

//_inlet_simple_types_tables_input_end
*/

int main()
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  // This is a generic basic SLIC logger
  axom::slic::SimpleLogger logger;

  std::string input = R""""(
    driver = {
        name = "Speed Racer",
        car = {
            make = "BestCompany",
            seats = 2,
            horsepower = 200
        }
    }
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
  // _inlet_simple_types_tables_add_start
  auto& driver_table = inlet.addTable("driver", "A description of driver");
  driver_table.addString("name", "Name of driver");

  auto& car_table = driver_table.addTable("car", "Car of driver");
  car_table.addString("make", "Make of car");
  car_table.addString("color", "Color of car").defaultValue("red");
  car_table.addInt("seats", "Number of seats");
  car_table.addInt("horsepower", "Amount of horsepower");
  // _inlet_simple_types_tables_add_end

  // Access values stored in the Datastore via Inlet
  // _inlet_simple_types_tables_access_start
  // Access values by fully qualified name from Inlet instance
  std::string name = inlet["driver/name"];

  // ... or... Get car table then access values from there
  auto car = inlet["driver/car"];
  std::string make = car["make"];
  std::string color = car["color"];
  int seats = car["seats"];
  int horsepower = car["horsepower"];
  // _inlet_simple_types_tables_access_end

  std::cout << "name = " << name << std::endl;
  std::cout << "make = " << make << std::endl;
  std::cout << "color = " << color << std::endl;
  std::cout << "seats = " << seats << std::endl;
  std::cout << "horsepower = " << horsepower << std::endl;

  return 0;
}
