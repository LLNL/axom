// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/io.hpp"

#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>

#include "conduit.hpp"

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
auto convertList(Node const &listToConvert, Converter converter)
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
std::vector<std::string> toStringList(Node const &listNode) {
    return convertList(listNode, std::mem_fn(&Node::as_string));
}

/**
 * Get the geometry specification for a shape.
 *
 * \param geometryNode the node describing the geometry
 * \return the geometry description for the shape
 */
Geometry getGeometry(Node const &geometryNode) {
    Geometry geometry;
    geometry.setFormat(geometryNode["format"].as_string());
    geometry.setPath(geometryNode["path"].as_string());
    return geometry;
}

/**
 * Convert a Node representing a Shape into a Shape object.
 *
 * \param shapeNode the shape as a conduit node
 * \return the shape as a Shape object
 */
Shape convertToShape(Node const &shapeNode) {
    Shape shape;
    shape.setName(shapeNode["name"].as_string());
    shape.setMaterial(shapeNode["material"].as_string());
    shape.setGeometry(getGeometry(shapeNode["geometry"]));

    if (shapeNode.has_child("replaces")) {
        if (shapeNode.has_child("does_not_replace")) {
            throw std::runtime_error("Can't have both 'replaces' and "
                                     "'does_not_replace' lists");
        }
        shape.setMaterialsReplaced(toStringList(shapeNode["replaces"]));
    } else if (shapeNode.has_child("does_not_replace")) {
        shape.setMaterialsNotReplaced(
                toStringList(shapeNode["does_not_replace"]));
    }

    return shape;
}

}

ShapeSet readShapeSet(std::istream &stream) {
    Node doc;
    std::string contents{std::istreambuf_iterator<char>(stream), {}};
    doc.parse(contents, "yaml");
    ShapeSet shapeSet;
    shapeSet.setShapes(convertList(doc["shapes"], convertToShape));
    return shapeSet;
}

}}
