// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"

#include <iostream>

#include "CLI11/CLI11.hpp"
#include "fmt/fmt.hpp"

#include "axom/slic/core/SimpleLogger.hpp"

namespace inlet = axom::inlet;
using Vector = inlet::FunctionType::Vector;

// Used to convert a std::vector (the type Inlet uses to represent arrays)
// to a geometric vector type
Vector toVector(const std::vector<double>& vec)
{
  // Narrow from std::size_t to int
  const int size = vec.size();
  return {{vec.data(), size}, size};
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
  Vector origin;
  enum class Type
  {
    Translate,
    Rotate,
    Slice
  } type;

  // Again, the union of the necessary members are defined as part of the schema
  static void defineSchema(inlet::Table& table)
  {
    // The rotation is in degrees
    table.addDouble("rotate").range(-180, 180);
    // Vectors are defined as arrays of doubles
    table.addDoubleArray("axis");
    table.addDoubleArray("center");
    table.addDoubleArray("translate");
    // The slice operation can have sub-entries, so we represent it as a struct
    // Note that Inlet does not require a 1-1 correspondence (or any correspondence)
    // between structures defined in the schema via addStruct and structures extracted
    // from the input file with FromInlet specializations (defined below)
    auto& slice = table.addStruct("slice");
    slice.addDouble("x");
    slice.addDouble("y");
    slice.addDouble("z");
    slice.addDoubleArray("origin");

    // Verify that exactly one type of operator is defined
    table.registerVerifier([](const inlet::Table& table) {
      const bool is_translate = table.contains("translate");
      const bool is_rotate = table.contains("rotate");
      const bool is_slice = table.contains("slice");

      // There can be only one
      if((is_translate && is_rotate) || (is_translate && is_slice) ||
         (is_rotate && is_slice))
      {
        return false;
      }

      return is_translate || is_rotate || is_slice;
    });
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
      result.translate = toVector(base["translate"]);
    }
    else if(base.contains("rotate"))
    {
      result.type = Operator::Type::Rotate;
      result.rotate = base["rotate"];
      result.axis = toVector(base["axis"]);
      result.center = toVector(base["center"]);
    }
    else if(base.contains("slice"))
    {
      result.type = Operator::Type::Slice;
      // Grab the substructure corresponding to the slice operation
      // and use it to populate the Operator instance
      auto slice = base["slice"];
      if(slice.contains("x"))
      {
        result.x = slice["x"];
      }
      if(slice.contains("y"))
      {
        result.y = slice["y"];
      }
      if(slice.contains("z"))
      {
        result.z = slice["z"];
      }
      if(slice.contains("origin"))
      {
        result.origin = toVector(slice["origin"]);
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
    os << fmt::format("       with origin: {0}\n", op.origin);
    break;
  default:
    SLIC_ERROR("Operator had unknown type");
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
    table.addString("units").defaultValue("cm").validValues({"cm", "m"});
    table.addInt("start_dimensions").defaultValue(3);
    // addStructArray is used to represent std::vector<T> where
    // T is any non-primitive type
    auto& ops_schema = table.addStructArray("operators");
    Operator::defineSchema(ops_schema);
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

    std::string units = base["units"];
    if(units == "cm")
    {
      result.units = Geometry::Units::Centimeters;
    }
    else if(units == "m")
    {
      result.units = Geometry::Units::Meters;
    }

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
    table.addString("name").required();
    table.addString("material").validValues({"steel", "wood", "plastic"});
    // addStruct is used for a single instance of a user-defined type
    auto& geom_schema = table.addStruct("geometry");
    Geometry::defineSchema(geom_schema);
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
    else if(material == "wood")
    {
      result.material = Shape::Material::Wood;
    }
    else if(material == "plastic")
    {
      result.material = Shape::Material::Plastic;
    }
    result.geom = base["geometry"].get<Geometry>();
    return result;
  }
};

// "rubber" is not one of the allowed shape materials, so this will fail verification
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
            origin: [40, 50, 60]
        - translate: [5, 6]
  - name: tire
    material: rubber
    geometry:
      format: test_format
      path: path/to/tire.format
      units: cm
      operators:
        - translate: [10, 20]
)";

int main(int argc, char** argv)
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  axom::slic::SimpleLogger logger;

  CLI::App app {"Example of Axom's Inlet component for nested structures"};
  bool docsEnabled {false};
  app.add_flag("--docs", docsEnabled, "Enables documentation generation");
  CLI11_PARSE(app, argc, argv);

  axom::sidre::DataStore ds;
  auto reader = std::unique_ptr<inlet::YAMLReader>(new inlet::YAMLReader());
  reader->parseString(input);
  inlet::Inlet inlet(std::move(reader), ds.getRoot());

  auto& shapes_table = inlet.addStructArray("shapes");
  Shape::defineSchema(shapes_table);

  if(inlet.verify())
  {
    SLIC_INFO("Verification was successful");
  }
  else
  {
    SLIC_INFO("Verification was unsuccessful");
  }

  auto shapes = inlet["shapes"].get<std::unordered_map<int, Shape>>();
  for(const auto& entry : shapes)
  {
    std::cout << entry.second << "\n";
  }

  if(docsEnabled)
  {
    std::unique_ptr<inlet::SphinxWriter> writer(
      new inlet::SphinxWriter("nested_structs.rst"));
    inlet.registerWriter(std::move(writer));
    inlet.writeDoc();
  }
}
