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
#include "axom/klee/IOUtil.hpp"

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
 * \param namedOperators any named operators that were parsed from the file
 * \return the geometry description for the shape
 */
Geometry getGeometry(const Node &geometryNode, Dimensions initialDimensions,
        const internal::NamedOperatorMap &namedOperators) {
    Geometry geometry;
    geometry.setInitialDimensions(initialDimensions);
    geometry.setFormat(geometryNode["format"].as_string());
    geometry.setPath(geometryNode["path"].as_string());
    if (geometryNode.has_child("initial_dimensions")) {
        geometry.setInitialDimensions(internal::toDimensions(
                geometryNode["initial_dimensions"]));
    }

    if (geometryNode.has_child("operators")) {
        auto operators = internal::parseGeometryOperators(
                geometryNode["operators"], geometry.getInitialDimensions(),
                namedOperators);
        geometry.setGeometryOperator(operators);
    }
    return geometry;
}

/**
 * Convert a Node representing a Shape into a Shape object.
 *
 * \param shapeNode the shape as a conduit node
 * \param namedOperators any named operators that were parsed from the file
 * \return the shape as a Shape object
 */
Shape convertToShape(const Node &shapeNode, Dimensions fileDimensions,
        const internal::NamedOperatorMap &namedOperators) {
    Shape shape;
    shape.setName(shapeNode["name"].as_string());
    shape.setMaterial(shapeNode["material"].as_string());
    shape.setGeometry(getGeometry(shapeNode["geometry"], fileDimensions,
            namedOperators));

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
 * Get all named geometry operators from the file
 * \param doc the node for the top-level document
 * \param initialDimensions the number of dimensions that operators should
 * start at unless otherwise specified
 * \return all named operators read from the document
 */
internal::NamedOperatorMap getNamedOperators(const Node &doc,
        Dimensions initialDimensions) {
    if (doc.has_child("named_operators")) {
        return internal::parseNamedGeometryOperators(doc["named_operators"],
                initialDimensions);
    }
    return {};
}
}

ShapeSet readShapeSet(std::istream &stream) {
    Node doc;
    std::string contents{std::istreambuf_iterator<char>(stream), {}};
    doc.parse(contents, "yaml");
    ShapeSet shapeSet;
    Dimensions dimensions = internal::toDimensions(doc["dimensions"]);
    auto namedOperators = getNamedOperators(doc, dimensions);
    shapeSet.setShapes(convertList(doc["shapes"],
            [dimensions, &namedOperators] (const Node &shapeNode) -> Shape {
                return convertToShape(shapeNode, dimensions, namedOperators);
            }));
    return shapeSet;
}

}}
