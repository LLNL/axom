// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"

#include <iostream>
#include <unordered_map>
#include "axom/slic/core/SimpleLogger.hpp"

using axom::inlet::FunctionType;
using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;

namespace inlet = axom::inlet;
// _inlet_userdef_simple_start
struct Car
{
  std::string make;
  std::string color;
  int seats;
  int horsepower;

  // A function can be used to define a schema for a particular struct
  // For convenience, it can be implemented as a static method of the struct
  static void defineSchema(inlet::Table& car_schema)
  {
    car_schema.addString("make", "Make of car");
    car_schema.addString("color", "Color of car").defaultValue("red");
    car_schema.addInt("seats", "Number of seats").range(2, 7);
    car_schema.addInt("horsepower", "Amount of horsepower");
  }
};
// _inlet_userdef_simple_end

// Additionally, each class should specialize this struct as follows
// in the global namespace so that Inlet can access it
// _inlet_userdef_simple_frominlet_start
template <>
struct FromInlet<Car>
{
  Car operator()(const inlet::Table& input_data)
  {
    Car result;
    result.make = input_data["make"];
    result.color = input_data["color"];
    result.seats = input_data["seats"];
    result.horsepower = input_data["horsepower"];
    return result;
  }
};
// _inlet_userdef_simple_frominlet_end

const std::string input = R"(
  car = {
      make = "BestCompany",
      seats = 2,
      horsepower = 200
  }
)";

int main()
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  axom::slic::SimpleLogger logger;

  DataStore ds;
  auto lr = std::make_unique<LuaReader>();
  lr->parseString(input);
  Inlet inlet(std::move(lr), ds.getRoot());

  // Create a table off the global table for the car object
  // then define its schema
  auto& car_schema = inlet.addStruct("car", "Vehicle description");
  Car::defineSchema(car_schema);

  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  // Extract the car object
  Car car = inlet["car"].get<Car>();
  std::cout << "make = " << car.make << std::endl;
  std::cout << "color = " << car.color << std::endl;
  std::cout << "seats = " << car.seats << std::endl;
  std::cout << "horsepower = " << car.horsepower << std::endl;
}
