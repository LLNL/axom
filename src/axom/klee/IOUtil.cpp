// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/IOUtil.hpp"

#include <stdexcept>

namespace axom
{
namespace klee
{
namespace internal
{
std::vector<double> toDoubleVector(const conduit::Node &listNode,
                                   std::size_t expectedSize)
{
  if(!listNode.dtype().is_number())
  {
    std::ostringstream message;
    message << "Not an array of numbers: " << listNode.path();
    throw std::invalid_argument(message.str());
  }

  conduit::Node valueAsDoubleArray;
  listNode.to_double_array(valueAsDoubleArray);
  double *values = valueAsDoubleArray.as_double_ptr();
  std::vector<double> result(valueAsDoubleArray.dtype().number_of_elements());
  std::copy(values, values + result.size(), result.begin());

  if(result.size() != expectedSize)
  {
    std::ostringstream message;
    message << listNode.name() << " should be a list of " << expectedSize
            << " numbers";
    throw std::invalid_argument(message.str());
  }
  return result;
}

double toDouble(const conduit::Node &value)
{
  auto &type = value.dtype();
  if(!type.is_number() || type.number_of_elements() != 1)
  {
    std::ostringstream message;
    message << value.name() << " should be a single number";
    throw std::invalid_argument(message.str());
  }
  return value.to_double();
}

Dimensions toDimensions(const conduit::Node &dimensionsNode)
{
  int dimensions = dimensionsNode.to_int();
  if(dimensions == 2)
  {
    return Dimensions::Two;
  }
  else if(dimensions == 3)
  {
    return Dimensions::Three;
  }
  throw std::invalid_argument("'dimensions' must be either 2 or 3");
}

std::tuple<LengthUnit, LengthUnit> getOptionalStartAndEndUnits(
  const conduit::Node &node)
{
  bool hasStartUnits = node.has_child("start_units");
  bool hasEndUnits = node.has_child("end_units");
  if(node.has_child("units"))
  {
    if(hasStartUnits || hasEndUnits)
    {
      throw std::invalid_argument(
        "Can't specify 'units' with 'start_units' or 'end_units'");
    }
    auto units = parseLengthUnits(node["units"].as_string());
    return std::make_tuple(units, units);
  }
  else if(hasStartUnits || hasEndUnits)
  {
    if(!(hasStartUnits && hasEndUnits))
    {
      throw std::invalid_argument(
        "Must specify both 'start_units' and 'end_units'");
    }
    auto startUnits = parseLengthUnits(node["start_units"].as_string());
    auto endUnits = parseLengthUnits(node["end_units"].as_string());
    return std::make_tuple(startUnits, endUnits);
  }
  return std::make_tuple(LengthUnit::unspecified, LengthUnit::unspecified);
}

std::tuple<LengthUnit, LengthUnit> getStartAndEndUnits(const conduit::Node &node)
{
  auto units = getOptionalStartAndEndUnits(node);
  if(std::get<0>(units) == LengthUnit::unspecified)
  {
    throw std::invalid_argument("Did not specify units");
  }
  return units;
}

}  // namespace internal
}  // namespace klee
}  // namespace axom