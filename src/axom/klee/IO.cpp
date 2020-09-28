// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/IO.hpp"

#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>

#include "conduit.hpp"

#include "axom/klee/GeometryOperatorsIO.hpp"

namespace axom { namespace klee {

namespace {

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
template<typename Converter>
auto convertList(const Node &listToConvert, Converter converter)
-> std::vector<decltype(converter(listToConvert.child(0)))> {
    std::vector<decltype(converter(listToConvert.child(0)))> converted;
    auto childIter = listToConvert.children();
    while (childIter.has_next()) {
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
std::vector<std::string> toStringList(const Node &listNode) {
    return convertList(listNode, std::mem_fn(&Node::as_string));
}

/**
 * Get the geometry specification for a shape.
 *
 * \param geometryNode the node describing the geometry
 * \param initialDimensions the initial dimensions of the shape
 * \return the geometry description for the shape
 */
Geometry getGeometry(const Node &geometryNode, int initialDimensions) {
    Geometry geometry;
    geometry.setInitialDimensions(initialDimensions);
    geometry.setFormat(geometryNode["format"].as_string());
    geometry.setPath(geometryNode["path"].as_string());
    if (geometryNode.has_child("initial_dimensions")) {
        geometry.setInitialDimensions(
                geometryNode["initial_dimensions"].to_int());
    }
    if (geometryNode.has_child("operators")) {
        auto operators = parseGeometryOperators(geometryNode["operators"],
                geometry.getInitialDimensions());
        geometry.setGeometryOperator(operators);
    }
    return geometry;
}

/**
 * Convert a Node representing a Shape into a Shape object.
 *
 * \param shapeNode the shape as a conduit node
 * \return the shape as a Shape object
 */
Shape convertToShape(const Node &shapeNode,
        int fileDimensions) {
    Shape shape;
    shape.setName(shapeNode["name"].as_string());
    shape.setMaterial(shapeNode["material"].as_string());
    shape.setGeometry(getGeometry(shapeNode["geometry"], fileDimensions));

    if (shapeNode.has_child("replaces")) {
        if (shapeNode.has_child("does_not_replace")) {
            throw std::invalid_argument("Can't have both 'replaces' and "
                                        "'does_not_replace' lists");
        }
        shape.setMaterialsReplaced(toStringList(shapeNode["replaces"]));
    } else if (shapeNode.has_child("does_not_replace")) {
        shape.setMaterialsNotReplaced(
                toStringList(shapeNode["does_not_replace"]));
    }

    return shape;
}

/**
 * Get the number of dimensions of shapes specified in the document.
 *
 * \param doc the document
 * \return the number of dimensions. Always 2 or 3.
 * \throws std::invalid_argument if the number of dimensions is invalid
 */
int getDimensions(const Node &doc) {
    int dimensions = doc["dimensions"].to_int();
    if (dimensions != 2 && dimensions != 3) {
        throw std::invalid_argument("'dimensions' must be either 2 or 3");
    }
    return dimensions;
}

}

ShapeSet readShapeSet(std::istream &stream) {
    Node doc;
    std::string contents{std::istreambuf_iterator<char>(stream), {}};
    doc.parse(contents, "yaml");
    ShapeSet shapeSet;
    int dimensions = getDimensions(doc);
    shapeSet.setShapes(convertList(doc["shapes"],
            [dimensions] (const Node &shapeNode) -> Shape {
                return convertToShape(shapeNode, dimensions);
            }));
    return shapeSet;
}

}}
