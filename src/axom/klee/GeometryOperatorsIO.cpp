// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/GeometryOperatorsIO.hpp"

#include "axom/klee/Geometry.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/IOUtil.hpp"

#include <functional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

namespace axom
{
namespace klee
{
namespace internal
{
namespace
{
using OpPtr = CompositeOperator::OpPtr;
using OperatorParser =
  std::function<OpPtr(const conduit::Node &, const TransformableGeometryProperties &)>;
using internal::toDouble;
using internal::toDoubleVector;
using primal::Point3D;
using primal::Vector3D;
using FieldSet = std::unordered_set<std::string>;

/**
 * Convert a Dimensions value to an integer.
 *
 * \param dims the number of dimensions
 * \return the number as an integer
 */
int toInt(Dimensions dims)
{
  if(dims == Dimensions::Two)
  {
    return 2;
  }
  return 3;
}

/**
 * Verify that a node has the correct fields.
 *
 * \param node the node to test
 * \param name the name of the node. This must be one of its fields
 * \param additionalRequiredFields any additional required fields
 * \param optionalFields any additional optional fields
 */
void verifyObjectFields(const conduit::Node &node,
                        const std::string &name,
                        const FieldSet &additionalRequiredFields,
                        const FieldSet &optionalFields)
{
  std::unordered_set<std::string> requiredParameters {additionalRequiredFields};
  requiredParameters.insert(name);

  for(auto &requiredParameter : requiredParameters)
  {
    if(!node.has_child(requiredParameter))
    {
      std::string message = "Missing required parameter \"";
      message += requiredParameter;
      message += "\" for operator \"";
      message += name;
      message += '"';
      throw std::invalid_argument(message);
    }
  }

  for(auto &parameter : node.child_names())
  {
    if(requiredParameters.count(parameter) == 0 &&
       optionalFields.count(parameter) == 0)
    {
      std::string message = "Unknown parameter \"";
      message += parameter;
      message += "\" for operator \"";
      message += name;
      message += '"';
      throw std::invalid_argument(message);
    }
  }
}

/**
 * Make an array-like structure from the given node.
 *
 * \tparam T the type of the array-like structure. Should be either
 * Point3D or Vector3D.
 * \param node the node to convert
 * \param dims the expected number of dimensions
 * \return the created value, scaled based on this parser's units
 */
template <typename T>
T toArrayType(const conduit::Node &node, Dimensions dims)
{
  static_assert(
    std::is_same<T, Point3D>::value || std::is_same<T, Vector3D>::value,
    "Type must be a Point3D or Vector3D");
  auto values = toDoubleVector(node, static_cast<std::size_t>(toInt(dims)));
  return T {values.data(), toInt(dims)};
}

/**
 * Make an array-like structure from the specified child of the
 * given node, if it exists.
 *
 * \tparam T the type of the array-like structure. Should be either
 * Point3D or Vector3D.
 * \param parent the parent node
 * \param name the name of the child
 * \param dims the expected number of dimensions
 * \return the created value, scaled based on this parser's units, or
 * the default value (not scaled)
 */
template <typename T>
T toOptionalArrayType(const conduit::Node &parent,
                      const std::string &name,
                      Dimensions dims,
                      T defaultValue)
{
  if(!parent.has_child(name))
  {
    return defaultValue;
  }
  return toArrayType<T>(parent[name], dims);
}

/**
 * Parse a "translate" operator.
 *
 * \param node the node from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseTranslate(const conduit::Node &node,
                     const TransformableGeometryProperties &startProperties)
{
  verifyObjectFields(node, "translate", FieldSet {}, FieldSet {});
  return std::make_shared<Translation>(
    toArrayType<Vector3D>(node["translate"], startProperties.dimensions),
    startProperties);
}

/**
 * Parse a "rotate" operator.
 *
 * \param node the node from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseRotate(const conduit::Node &node,
                  const TransformableGeometryProperties &startProperties)
{
  if(startProperties.dimensions == Dimensions::Two)
  {
    verifyObjectFields(node, "rotate", FieldSet {}, {"center"});
    Vector3D axis {0, 0, 1};
    return std::make_shared<Rotation>(
      toDouble(node["rotate"]),
      toOptionalArrayType<Point3D>(node,
                                   "center",
                                   startProperties.dimensions,
                                   {0, 0, 0}),
      axis,
      startProperties);
  }
  else
  {
    verifyObjectFields(node, "rotate", {"axis"}, {"center"});
    return std::make_shared<Rotation>(
      toDouble(node["rotate"]),
      toOptionalArrayType<Point3D>(node,
                                   "center",
                                   startProperties.dimensions,
                                   {0, 0, 0}),
      toArrayType<Vector3D>(node["axis"], Dimensions::Three),
      startProperties);
  }
}

/**
 * Make a slice, ensuring all the values are valid.
 *
 * \param origin the origin of the coordinate system
 * \param normal a vector normal to the plane
 * \param up a vector which defines the positive Y direction
 * \param startProperties the properties before the slice
 * \return the created operator
 * \throws std::invalid_argument if any value is invalid
 */
OpPtr makeCheckedSlice(Point3D origin,
                       Vector3D normal,
                       Vector3D up,
                       const TransformableGeometryProperties &startProperties)
{
  if(normal.is_zero())
  {
    throw std::invalid_argument(
      "The 'normal' vector must not be a zero vector");
  }
  if(!utilities::isNearlyEqual(normal.dot(up), 0.0))
  {
    throw std::invalid_argument(
      "The 'normal' and 'up' vectors must be perpendicular");
  }
  return std::make_shared<SliceOperator>(origin, normal, up, startProperties);
}

/**
 * Get the origin to use for a perpendicular slice.
 *
 * \param sliceNode the node describing the slice
 * \param planeName the name of the plane ("x", "y", or "z")
 * \param defaultNormal the default normal vector
 * \return the point to use as the origin
 */
primal::Point3D getPerpendicularSliceOrigin(const conduit::Node &sliceNode,
                                            char const *planeName,
                                            const primal::Vector3D &defaultNormal)
{
  double axisIntercept = toDouble(sliceNode[planeName]);

  primal::Point3D defaultOrigin;
  int nonZeroIndex = -1;
  for(int i = 0; i < 3; ++i)
  {
    defaultOrigin[i] = axisIntercept * defaultNormal[i];
    if(!utilities::isNearlyEqual(defaultNormal[i], 0.0))
    {
      nonZeroIndex = i;
    }
  }

  if(!sliceNode.has_child("origin"))
  {
    return defaultOrigin;
  }

  primal::Point3D givenOrigin =
    toArrayType<Point3D>(sliceNode["origin"], Dimensions::Three);
  if(givenOrigin[nonZeroIndex] != axisIntercept)
  {
    throw std::invalid_argument("The origin must be on the slice plane");
  }
  return givenOrigin;
}

/**
 * Get the normal vector to use for a perpendicular slice.
 *
 * \param sliceNode the node describing the slice
 * \param defaultNormal the default normal vector
 * \return the vector to use as the normal
 */
primal::Vector3D getPerpendicularSliceNormal(const conduit::Node &sliceNode,
                                             const primal::Vector3D &defaultNormal)
{
  if(!sliceNode.has_child("normal"))
  {
    return defaultNormal;
  }

  primal::Vector3D givenNormal =
    toArrayType<Vector3D>(sliceNode["normal"], Dimensions::Three);
  auto cross = primal::Vector3D::cross_product(givenNormal, defaultNormal);
  bool parallel = cross.is_zero();
  if(!parallel)
  {
    throw std::invalid_argument("Invalid normal");
  }
  return givenNormal;
}

/**
 * Read a perpendicular slice.
 *
 * \param sliceNode the node describing the slice
 * \param planeName the name of the plane ("x", "y", or "z")
 * \param defaultNormal the default normal vector for the type of plane
 *  being parsed
 * \param defaultUp the default up vector for the plane being parsed
 * \param startProperties the properties prior to this operator
 * \return the parsed plane
 */
OpPtr readPerpendicularSlice(const conduit::Node &sliceNode,
                             char const *planeName,
                             Vector3D const &defaultNormal,
                             Vector3D const &defaultUp,
                             const TransformableGeometryProperties &startProperties)
{
  verifyObjectFields(sliceNode,
                     planeName,
                     FieldSet {},
                     {"origin", "normal", "up"});
  const primal::Vector3D defaultNormalVec {defaultNormal.data()};

  auto origin =
    getPerpendicularSliceOrigin(sliceNode, planeName, defaultNormalVec);
  auto normal = getPerpendicularSliceNormal(sliceNode, defaultNormalVec);
  auto up =
    toOptionalArrayType<Vector3D>(sliceNode, "up", Dimensions::Three, defaultUp);

  return makeCheckedSlice(origin, normal, up, startProperties);
}

/**
 * Parse a "slice" operator.
 *
 * \param node the node from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseSlice(const conduit::Node &node,
                 const TransformableGeometryProperties &startProperties)
{
  if(startProperties.dimensions != Dimensions::Three)
  {
    throw std::invalid_argument("Cannot do a slice from 2D");
  }
  verifyObjectFields(node, "slice", FieldSet {}, FieldSet {});
  auto &properties = node["slice"];
  if(properties.has_child("x"))
  {
    return readPerpendicularSlice(properties,
                                  "x",
                                  {1, 0, 0},
                                  {0, 0, 1},
                                  startProperties);
  }
  else if(properties.has_child("y"))
  {
    return readPerpendicularSlice(properties,
                                  "y",
                                  {0, 1, 0},
                                  {1, 0, 0},
                                  startProperties);
  }
  else if(properties.has_child("z"))
  {
    return readPerpendicularSlice(properties,
                                  "z",
                                  {0, 0, 1},
                                  {0, 1, 0},
                                  startProperties);
  }

  verifyObjectFields(properties, "origin", {"normal", "up"}, FieldSet {});
  return makeCheckedSlice(
    toArrayType<Point3D>(properties["origin"], Dimensions::Three),
    toArrayType<Vector3D>(properties["normal"], Dimensions::Three),
    toArrayType<Vector3D>(properties["up"], Dimensions::Three),
    startProperties);
}

/**
 * Parse a "scale" operator.
 *
 * \param node the node from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseScale(const conduit::Node &node,
                 const TransformableGeometryProperties &startProperties)
{
  verifyObjectFields(node, "scale", FieldSet {}, FieldSet {});
  auto &scaleNode = node["scale"];
  if(scaleNode.dtype().is_number() && scaleNode.dtype().number_of_elements() == 1)
  {
    double factor = scaleNode.to_double();
    return std::make_shared<Scale>(factor, factor, factor, startProperties);
  }
  auto factors = toDoubleVector(scaleNode, toInt(startProperties.dimensions));
  if(startProperties.dimensions == Dimensions::Two)
  {
    factors.emplace_back(1.0);
  }
  return std::make_shared<Scale>(factors[0],
                                 factors[1],
                                 factors[2],
                                 startProperties);
}

/**
 * Parse a "convert_units_to" operator.
 *
 * \param node the node from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseConvertUnits(const conduit::Node &node,
                        const TransformableGeometryProperties &startProperties)
{
  verifyObjectFields(node, "convert_units_to", FieldSet {}, FieldSet {});
  auto endUnits = parseLengthUnits(node["convert_units_to"].as_string());
  return std::make_shared<UnitConverter>(endUnits, startProperties);
}

/**
 * Parse an operator specified via the "ref" command.
 *
 * \param node the node from which to read the operator
 * \param startProperties the properties before the "ref" command
 * \param namedOperators a map of named operators from which to get
 * referenced operators
 * \return the created operator
 */
OpPtr parseRef(const conduit::Node &node,
               const TransformableGeometryProperties &startProperties,
               const NamedOperatorMap &namedOperators)
{
  verifyObjectFields(node, "ref", FieldSet {}, FieldSet {});
  std::string const &operatorName = node["ref"].as_string();
  auto opIter = namedOperators.find(operatorName);
  if(opIter == namedOperators.end())
  {
    std::string message = "No operator named '";
    message += operatorName;
    message += '\'';
    throw std::invalid_argument(message);
  }
  auto referencedOperator = opIter->second;
  bool startUnitsMatch =
    startProperties.units == referencedOperator->getStartProperties().units;
  bool endUnitsMatch =
    startProperties.units == referencedOperator->getEndProperties().units;

  if(startUnitsMatch && endUnitsMatch)
  {
    return referencedOperator;
  }

  auto compositeWithConversions =
    std::make_shared<CompositeOperator>(startProperties);
  if(!startUnitsMatch)
  {
    compositeWithConversions->addOperator(std::make_shared<UnitConverter>(
      referencedOperator->getStartProperties().units,
      startProperties));
  }
  compositeWithConversions->addOperator(referencedOperator);
  if(!endUnitsMatch)
  {
    compositeWithConversions->addOperator(
      std::make_shared<UnitConverter>(startProperties.units,
                                      referencedOperator->getEndProperties()));
  }
  return compositeWithConversions;
}

/**
 * Parse a single operator
 *
 * \param node the node from which to parse the operator
 * \param startProperties the properties before the operator
 * \param namedOperators a map of named operators from which to get
 * referenced operators
 * \return the created operator
 */
OpPtr parseSingleOperator(const conduit::Node &node,
                          const TransformableGeometryProperties &startProperties,
                          const NamedOperatorMap &namedOperators)
{
  std::unordered_map<std::string, OperatorParser> parsers {
    {"translate", parseTranslate},
    {"rotate", parseRotate},
    {"slice", parseSlice},
    {"scale", parseScale},
    {"convert_units_to", parseConvertUnits},
    {"ref",
     [&namedOperators](const conduit::Node &opNode,
                       const TransformableGeometryProperties &startProperties) {
       return parseRef(opNode, startProperties, namedOperators);
     }},
  };

  for(auto &entry : parsers)
  {
    if(node.has_child(entry.first))
    {
      return entry.second(node, startProperties);
    }
  }

  std::string message = "Invalid transformation: \n";
  message += node.to_json_default();
  throw std::invalid_argument(message);
}

}  // namespace

std::shared_ptr<const GeometryOperator> parseGeometryOperators(
  const conduit::Node &node,
  const TransformableGeometryProperties &startProperties,
  const NamedOperatorMap &namedOperators)
{
  auto composite = std::make_shared<CompositeOperator>(startProperties);
  TransformableGeometryProperties currentProperties = startProperties;
  for(auto iter = node.children(); iter.has_next();)
  {
    auto op = parseSingleOperator(iter.next(), currentProperties, namedOperators);
    currentProperties = op->getEndProperties();
    composite->addOperator(op);
  }
  return composite;
}

NamedOperatorMap parseNamedGeometryOperators(const conduit::Node &node,
                                             Dimensions startDimensions)
{
  NamedOperatorMap namedOperators;

  for(auto iter = node.children(); iter.has_next();)
  {
    auto &namedOperator = iter.next();
    std::string name = namedOperator["name"].as_string();

    Dimensions dimensions = startDimensions;
    if(namedOperator.has_child("start_dimensions"))
    {
      dimensions = toDimensions(namedOperator["start_dimensions"]);
    }

    auto units = internal::getStartAndEndUnits(namedOperator);

    TransformableGeometryProperties startProperties {
      dimensions,
      std::get<0>(units),
    };
    auto op = parseGeometryOperators(namedOperator["value"],
                                     startProperties,
                                     namedOperators);

    if(op->getEndProperties().units != std::get<1>(units))
    {
      throw std::invalid_argument(
        "Specified end units did not match actual units");
    }
    namedOperators.insert({name, op});
  }
  return namedOperators;
}

}  // namespace internal
}  // namespace klee
}  // namespace axom
