// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/IO.hpp"

#include <fstream>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>

#include "conduit.hpp"

#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/GeometryOperatorsIO.hpp"
#include "axom/klee/IOUtil.hpp"

namespace axom
{
namespace klee
{
namespace
{
using conduit::Node;

/**
 * Convert a conduit Node that represents a list to a vector
 *
 * \tparam Converter a type of a function that can convert a child node
 * to the desired type
 * \param listToConvert the node that represents the list to convert
 * \param converter a function which can convert a child node to the desired
 * value
 * \return a vector of converted values
 */
template <typename Converter>
auto convertList(const Node &listToConvert, Converter converter)
  -> std::vector<decltype(converter(listToConvert.child(0)))>
{
  std::vector<decltype(converter(listToConvert.child(0)))> converted;
  auto childIter = listToConvert.children();
  while(childIter.has_next())
  {
    converted.emplace_back(converter(childIter.next()));
  }
  return converted;
}

/**
 * Convert a Node which is supposed to be an array of strings to the equivalent
 * vector.
 *
 * \param listNode the node that represents the list
 * \return  the node as a string vector
 */
std::vector<std::string> toStringList(const Node &listNode)
{
  return convertList(listNode, std::mem_fn(&Node::as_string));
}

/**
 * Get the geometry specification for a shape.
 *
 * \param geometryNode the node describing the geometry
 * \param fileDimensions the dimensions of the file
 * \param namedOperators any named operators that were parsed from the file
 * \return the geometry description for the shape
 */
Geometry getGeometry(const Node &geometryNode,
                     Dimensions fileDimensions,
                     const internal::NamedOperatorMap &namedOperators)
{
  Dimensions startDimensions = fileDimensions;
  if(geometryNode.has_child("start_dimensions"))
  {
    startDimensions =
      internal::toDimensions(geometryNode["start_dimensions"]);
  }

  LengthUnit startUnits, endUnits;
  std::tie(startUnits, endUnits) =
    internal::getOptionalStartAndEndUnits(geometryNode);
  TransformableGeometryProperties startProperties {
    startDimensions,
    startUnits,
  };
  TransformableGeometryProperties endProperties = startProperties;

  std::shared_ptr<GeometryOperator const> operator_;
  if(geometryNode.has_child("operators"))
  {
    if(startUnits == LengthUnit::unspecified)
    {
      throw std::invalid_argument(
        "Cannot specify operators without specifying units");
    }
    operator_ = internal::parseGeometryOperators(geometryNode["operators"],
                                                 startProperties,
                                                 namedOperators);
    endProperties = operator_->getEndProperties();
  }

  if(endProperties.dimensions != fileDimensions)
  {
    throw std::invalid_argument(
      "Did not end up in the number of dimensions specified by the file");
  }

  return Geometry {startProperties,
                   geometryNode["format"].as_string(),
                   geometryNode["path"].as_string(),
                   operator_};
}

/**
 * Convert a Node representing a Shape into a Shape object.
 *
 * \param shapeNode the shape as a conduit node
 * \param namedOperators any named operators that were parsed from the file
 * \return the shape as a Shape object
 */
Shape convertToShape(const Node &shapeNode,
                     Dimensions fileDimensions,
                     const internal::NamedOperatorMap &namedOperators)
{
  std::vector<std::string> materialsReplaced;
  std::vector<std::string> materialsNotReplaced;
  if(shapeNode.has_child("replaces"))
  {
    if(shapeNode.has_child("does_not_replace"))
    {
      throw std::invalid_argument(
        "Can't have both 'replaces' and "
        "'does_not_replace' lists");
    }
    materialsReplaced = toStringList(shapeNode["replaces"]);
  }
  else if(shapeNode.has_child("does_not_replace"))
  {
    materialsNotReplaced = toStringList(shapeNode["does_not_replace"]);
  }

  return Shape {
    shapeNode["name"].as_string(),
    shapeNode["material"].as_string(),
    materialsReplaced,
    materialsNotReplaced,
    getGeometry(shapeNode["geometry"], fileDimensions, namedOperators)};
}

/**
 * Get all named geometry operators from the file
 * \param doc the node for the top-level document
 * \param startDimensions the number of dimensions that operators should
 * start at unless otherwise specified
 * \return all named operators read from the document
 */
internal::NamedOperatorMap getNamedOperators(const Node &doc,
                                             Dimensions startDimensions)
{
  if(doc.has_child("named_operators"))
  {
    return internal::parseNamedGeometryOperators(doc["named_operators"],
                                                 startDimensions);
  }
  return internal::NamedOperatorMap {};
}
}  // namespace

ShapeSet readShapeSet(std::istream &stream)
{
  Node doc;
  std::string contents {std::istreambuf_iterator<char>(stream), {}};
  doc.parse(contents, "yaml");
  ShapeSet shapeSet;
  Dimensions dimensions = internal::toDimensions(doc["dimensions"]);
  auto namedOperators = getNamedOperators(doc, dimensions);
  shapeSet.setShapes(
    convertList(doc["shapes"],
                [=, &namedOperators](const Node &shapeNode) -> Shape {
                  return convertToShape(shapeNode, dimensions, namedOperators);
                }));
  return shapeSet;
}

ShapeSet readShapeSet(const std::string &filePath)
{
  std::ifstream fin {filePath};
  auto shapeSet = readShapeSet(fin);
  fin.close();
  shapeSet.setPath(filePath);
  return shapeSet;
}
}  // namespace klee
}  // namespace axom
