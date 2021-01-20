#include "axom/inlet.hpp"

#include "CLI11/CLI11.hpp"

#include "axom/slic/core/SimpleLogger.hpp"

namespace inlet = axom::inlet;
using Vector = inlet::FunctionType::Vec3D;

Vector mapToVector(const std::unordered_map<int, double>& map)
{
  Vector result;
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
  static void defineSchema(inlet::Table& table)
  {
    table.addDouble("rotate");
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
    else if(base.contains("slice"))
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
    auto ops = base["operators"].get<std::unordered_map<int, Operator>>();
    for(const auto& ele : ops)
    {
      result.operators.push_back(ele.second);
    }
    return result;
  }
};

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
    auto& geom_table = table.addTable("geometry");
    Geometry::defineSchema(geom_table);
  }
};

int main(int argc, char** argv)
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  axom::slic::SimpleLogger logger;

  CLI::App app {"Example of Axom's Inlet component with nested types"};
  std::string inputFileName;
  auto opt = app.add_option("--file", inputFileName, "Path to input file");
  opt->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  axom::sidre::DataStore ds;
  auto reader = std::make_unique<inlet::YAMLReader>();
  reader->parseFile(inputFileName);
  inlet::Inlet inlet(std::move(reader), ds.getRoot());
}
