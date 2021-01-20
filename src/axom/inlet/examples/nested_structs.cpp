// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"

#include <iostream>

#include "fmt/fmt.hpp"

#include "axom/slic/core/SimpleLogger.hpp"

namespace inlet = axom::inlet;
using Vector = inlet::FunctionType::Vec3D;

// Used to convert an unordered map (the type Inlet uses) to represent arrays
// to a geometric vector type
// FIXME: Clean this up/remove when PR #429 is merged
Vector mapToVector(const std::unordered_map<int, double>& map)
{
  Vector result;
  // Make sure the elements are accessed in order
  for(int i = 0; i < 3; i++)
  {
    const auto ele = map.find(i);
    if(ele != map.end())
    {
      result[i] = ele->second;
    }
  }
  return result;
}

// A union of the members required for each of the operations is stored for simplicity
struct Operator
{
  double rotate;
  Vector axis;
  Vector center;
  Vector translate;
  double x;
  double y;
  double z;
  enum class Type
  {
    Translate,
    Rotate,
    Slice
  } type;

  // Again, the union of the necessary members are defined as part of the schema
  static void defineSchema(inlet::Table& table)
  {
    table.addDouble("rotate");
    // Vectors are defined as arrays of doubles
    table.addDoubleArray("axis");
    table.addDoubleArray("center");
    table.addDoubleArray("translate");
    table.addDouble("x");
    table.addDouble("y");
    table.addDouble("z");
  }
};

template <>
struct FromInlet<Operator>
{
  Operator operator()(const inlet::Table& base)
  {
    Operator result;
    // Even though all the possible members are part of the schema, the
    // "contains" operation can be used to check what is actually present
    // in the input file
    if(base.contains("translate"))
    {
      result.type = Operator::Type::Translate;
      result.translate = mapToVector(base["translate"]);
    }
    else if(base.contains("rotate"))
    {
      result.type = Operator::Type::Rotate;
      result.rotate = base["rotate"];
      result.axis = mapToVector(base["axis"]);
      result.center = mapToVector(base["center"]);
    }
    else if(base.contains("x") || base.contains("y") || base.contains("z"))
    {
      result.type = Operator::Type::Slice;
      if(base.contains("x"))
      {
        result.x = base["x"];
      }
      if(base.contains("y"))
      {
        result.y = base["y"];
      }
      if(base.contains("z"))
      {
        result.z = base["z"];
      }
    }
    return result;
  }
};

std::ostream& operator<<(std::ostream& os, const Operator& op)
{
  switch(op.type)
  {
  case Operator::Type::Translate:
    os << "   Translation operator:\n";
    os << fmt::format("      with vector: {0}\n", op.translate);
    break;
  case Operator::Type::Rotate:
    os << "   Rotation operator:\n";
    os << fmt::format("      with axis: {0}\n", op.axis);
    os << fmt::format("      with center: {0}\n", op.center);
    break;
  case Operator::Type::Slice:
    os << "   Slice operator:\n";
    os << fmt::format("      with x-coord: {0}\n", op.x);
    os << fmt::format("      with y-coord: {0}\n", op.y);
    os << fmt::format("      with z-coord: {0}\n", op.z);
    break;
  }
  return os;
}

struct Geometry
{
  std::vector<Operator> operators;
  std::string format;
  std::string path;
  enum class Units
  {
    Centimeters,
    Meters
  } units;
  int start_dim;
  static void defineSchema(inlet::Table& table)
  {
    table.addString("format");
    table.addString("path");
    // A string is used to represent the enumeration
    table.addString("units");
    table.addInt("start_dimensions");
    // addGenericArray is used to represent std::vector<T> where
    // T is any non-primitive type
    auto& ops_table = table.addGenericArray("operators");
    Operator::defineSchema(ops_table);
  }
};

template <>
struct FromInlet<Geometry>
{
  Geometry operator()(const inlet::Table& base)
  {
    Geometry result;
    result.format = base["format"];
    result.path = base["path"];
    result.start_dim = base["start_dimensions"];

    // FIXME: Remove when PR #429 is merged
    auto ops = base["operators"].get<std::unordered_map<int, Operator>>();
    for(const auto& ele : ops)
    {
      result.operators.push_back(ele.second);
    }
    return result;
  }
};

std::ostream& operator<<(std::ostream& os, const Geometry& geom)
{
  os << fmt::format("Geometry in format: '{0}'\n", geom.format);
  os << fmt::format("  with path: '{0}'\n", geom.path);
  for(const auto& op : geom.operators)
  {
    os << op;
  }
  return os;
}

struct Shape
{
  std::string name;
  enum class Material
  {
    Steel,
    Wood,
    Plastic
  } material;
  Geometry geom;
  static void defineSchema(inlet::Table& table)
  {
    table.addString("name");
    table.addString("material");
    // addStruct is used for a single instance of a user-defined type
    auto& geom_table = table.addStruct("geometry");
    Geometry::defineSchema(geom_table);
  }
};

std::ostream& operator<<(std::ostream& os, const Shape& shape)
{
  os << fmt::format("Shape: '{0}'\n", shape.name);
  os << shape.geom;
  return os;
}

template <>
struct FromInlet<Shape>
{
  Shape operator()(const inlet::Table& base)
  {
    Shape result;
    result.name = base["name"];
    std::string material = base["material"];
    if(material == "steel")
    {
      result.material = Shape::Material::Steel;
    }
    result.geom = base["geometry"].get<Geometry>();
    return result;
  }
};

const std::string input = R"(
dimensions: 2
shapes:
  - name: wheel
    material: steel
    geometry:
      format: test_format
      path: path/to/file.format
      units: cm
      start_dimensions: 3
      operators:
        - translate: [10, 20, 30]
        - rotate: 45
          axis: [1, 2, 3]
          center: [4, 5, 6]
        - slice:
          x: 10
        - translate: [5, 6]
)";

int main()
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  axom::slic::SimpleLogger logger;

  axom::sidre::DataStore ds;
  auto reader = std::make_unique<inlet::YAMLReader>();
  reader->parseString(input);
  inlet::Inlet inlet(std::move(reader), ds.getRoot());

  auto& shapes_table = inlet.addGenericArray("shapes");
  Shape::defineSchema(shapes_table);

  auto shapes = inlet["shapes"].get<std::unordered_map<int, Shape>>();
  for(const auto& entry : shapes)
  {
    std::cout << entry.second << "\n";
  }
}
