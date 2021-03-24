// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/IOUtil.hpp"

#include "axom/inlet.hpp"
#include "axom/klee/KleeError.hpp"

namespace axom
{
namespace klee
{
namespace internal
{
std::vector<double> toDoubleVector(inlet::Proxy const &field,
                                   Dimensions expectedDims,
                                   char const *fieldName)
{
  auto expectedSize = static_cast<std::size_t>(expectedDims);
  auto values = field.get<std::vector<double>>();
  auto actualSize = values.size();
  if(actualSize != expectedSize)
  {
    throw KleeError(fmt::format("Wrong size for {}. Expected {}. Got {}.",
                                fieldName,
                                expectedSize,
                                actualSize));
  }
  return values;
}

template <typename T>
T toArrayLike(inlet::Proxy const &parent,
              char const *fieldName,
              Dimensions expectedDims)
{
  auto values = toDoubleVector(parent[fieldName], expectedDims, fieldName);
  return T {values.data(), static_cast<int>(expectedDims)};
}

template <typename T>
T toArrayLike(inlet::Proxy const &parent,
              char const *fieldName,
              Dimensions expectedDims,
              const T &defaultValue)
{
  if(parent.contains(fieldName))
  {
    return toArrayLike<T>(parent, fieldName, expectedDims);
  }
  return defaultValue;
}

primal::Point3D toPoint(inlet::Proxy const &parent,
                        char const *fieldName,
                        Dimensions expectedDims)
{
  return toArrayLike<primal::Point3D>(parent, fieldName, expectedDims);
}

primal::Point3D toPoint(inlet::Proxy const &parent,
                        char const *fieldName,
                        Dimensions expectedDims,
                        const primal::Point3D &defaultValue)
{
  return toArrayLike(parent, fieldName, expectedDims, defaultValue);
}

primal::Vector3D toVector(inlet::Proxy const &parent,
                          char const *fieldName,
                          Dimensions expectedDims)
{
  return toArrayLike<primal::Vector3D>(parent, fieldName, expectedDims);
}

primal::Vector3D toVector(inlet::Proxy const &parent,
                          char const *fieldName,
                          Dimensions expectedDims,
                          const primal::Vector3D &defaultValue)
{
  return toArrayLike(parent, fieldName, expectedDims, defaultValue);
}

std::tuple<LengthUnit, LengthUnit> getOptionalStartAndEndUnits(
  const inlet::Proxy &proxy)
{
  bool hasStartUnits = proxy.contains("start_units");
  bool hasEndUnits = proxy.contains("end_units");
  if(proxy.contains("units"))
  {
    if(hasStartUnits || hasEndUnits)
    {
      throw KleeError(
        "Can't specify 'units' with 'start_units' or 'end_units'");
    }
    auto units = parseLengthUnits(proxy["units"]);
    return std::make_tuple(units, units);
  }
  else if(hasStartUnits || hasEndUnits)
  {
    if(!(hasStartUnits && hasEndUnits))
    {
      throw KleeError("Must specify both 'start_units' and 'end_units'");
    }
    auto startUnits = parseLengthUnits(proxy["start_units"]);
    auto endUnits = parseLengthUnits(proxy["end_units"]);
    return std::make_tuple(startUnits, endUnits);
  }
  return std::make_tuple(LengthUnit::unspecified, LengthUnit::unspecified);
}

std::tuple<LengthUnit, LengthUnit> getStartAndEndUnits(const inlet::Proxy &proxy)
{
  auto units = getOptionalStartAndEndUnits(proxy);
  if(std::get<0>(units) == LengthUnit::unspecified)
  {
    throw KleeError("Did not specify units");
  }
  return units;
}

void defineUnitsSchema(inlet::Container &container,
                       const char *unitsDescription,
                       const char *startUnitsDescription,
                       const char *endUnitsDescription)
{
  container.addString("start_units", startUnitsDescription);
  container.addString("end_units", endUnitsDescription);
  container.addString("units", unitsDescription);
  // Don't do custom validator here because getOptionalStartAndEndUnits()
  // verifies the right combination is specified. If we were to add a
  // custom validator, we would have to repeat some of the logic when
  // figuring out which fields to use.
}

inlet::VerifiableScalar &defineDimensionsField(inlet::Container &parent,
                                               const char *name,
                                               const char *description)
{
  return parent.addInt(name, description).range(2, 3);
}

Dimensions toDimensions(const inlet::Proxy &dimProxy)
{
  return static_cast<Dimensions>(dimProxy.get<int>());
}

}  // namespace internal
}  // namespace klee
}  // namespace axom
