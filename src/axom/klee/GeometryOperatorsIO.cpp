// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/GeometryOperatorsIO.hpp"

#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/IOUtil.hpp"

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace axom { namespace klee { namespace internal {
namespace {

using OpPtr = CompositeOperator::OpPtr;
using OperatorParser = std::function<OpPtr(const conduit::Node &, Dimensions)>;
using internal::toDoubleVector;
using internal::toDouble;
using primal::Point3D;
using primal::Vector3D;
using FieldSet = std::unordered_set<std::string>;

/**
 * Convert a Dimensions value to an integer.
 *
 * \param dims the number of dimensions
 * \return the number as an integer
 */
int toInt(Dimensions dims) {
    if (dims == Dimensions::Two) {
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
        const FieldSet &optionalFields) {

    std::unordered_set<std::string> requiredParameters{
            additionalRequiredFields};
    requiredParameters.insert(name);

    for (auto &requiredParameter : requiredParameters) {
        if (!node.has_child(requiredParameter)) {
            std::string message = "Missing required parameter \"";
            message += requiredParameter;
            message += "\" for operator \"";
            message += name;
            message += '"';
            throw std::invalid_argument(message);
        }
    }

    for (auto &parameter : node.child_names()) {
        if (requiredParameters.count(parameter) == 0
            && optionalFields.count(parameter) == 0) {
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
 * Make an object of the specified type from a node that should be a list of
 * numbers.
 *
 * \tparam T the type of object to make. Should be an instantiation of
 * either primal::Vector<double, N> or primal::Point<double, N>.
 * \param node the Conduit node containing the object
 * \param dimensions the required number of dimensions. Does not have
 * to match N. Whatever numbers are found will be given to the object's
 * constructor.
 * \return the created object
 * \throws std::invalid_argument if the node is not an array with the
 * specified number of entries.
 */
template<typename T>
T makeFromDoubleVector(const conduit::Node &node, Dimensions dimensions) {
    auto values = toDoubleVector(node,
            static_cast<std::size_t>(toInt(dimensions)));
    return T{values.data(), toInt(dimensions)};
}

/**
 * Create a primal::Vector3D from the given node.
 *
 * \param node the Conduit node to convert
 * \param dimensions the expected number of entries in the vector. If less
 * than 3, trailing values will be zero.
 * \return the created vector
 */
Vector3D toVector3D(const conduit::Node &node, Dimensions dimensions) {
    return makeFromDoubleVector<Vector3D>(node, dimensions);
}

/**
 * Create a primal::Point3D from the given node.
 *
 * \param node the Conduit node to convert
 * \param dimensions the expected number of entries in the point. If less
 * than 3, trailing values will be zero.
 * \return the created point
 */
Point3D toPoint3D(const conduit::Node &node, Dimensions dimensions) {
    return makeFromDoubleVector<Point3D>(node, dimensions);
}

/**
 * Get an optional Point3D from a node.
 *
 * \param parent the parent node
 * \param name the name of the optional child node
 * \param dimensions the expected number of dimensions in the array
 * \param defaultValue the default value to use when the child is missing
 * \return the value of the child, or the default value if there is no such
 * child in the parent node.
 */
Point3D getOptionalPoint(const conduit::Node &parent, const std::string &name,
        Dimensions dimensions, const Point3D &defaultValue) {
    if (!parent.has_child(name)) {
        return defaultValue;
    }
    return toPoint3D(parent[name], dimensions);
}

/**
 * Get an optional Vector3D from a node.
 *
 * \param parent the parent node
 * \param name the name of the optional child node
 * \param dimensions the expected number of dimensions in the array
 * \param defaultValue the default value to use when the child is missing
 * \return the value of the child, or the default value if there is no such
 * child in the parent node.
 */
Vector3D getOptionalVector(const conduit::Node &parent, const std::string &name,
        Dimensions dimensions, const Vector3D &defaultValue) {
    if (!parent.has_child(name)) {
        return defaultValue;
    }
    return toPoint3D(parent[name], dimensions);
}

/**
 * Parse a "translate" operator.
 *
 * \param node the node from which to read the operator
 * \param dimensions the expected number of dimensions
 * \return the created operator
 */
OpPtr parseTranslate(const conduit::Node &node, Dimensions dimensions) {
    verifyObjectFields(node, "translate", FieldSet{}, FieldSet{});
    return std::make_shared<Translation>(
            toVector3D(node["translate"], dimensions),
            dimensions);
}

/**
 * Parse a "rotate" operator.
 *
 * \param node the node from which to read the operator
 * \param dimensions the expected number of dimensions
 * \return the created operator
 */
OpPtr parseRotate(const conduit::Node &node, Dimensions dimensions) {
    if (dimensions == Dimensions::Two) {
        verifyObjectFields(node, "rotate", FieldSet{}, {"center"});
        Vector3D axis{0, 0, 1};
        return std::make_shared<Rotation>(
                toDouble(node["rotate"]),
                getOptionalPoint(node, "center", dimensions, {0, 0, 0}),
                axis, dimensions);
    } else {
        verifyObjectFields(node, "rotate", {"axis"}, {"center"});
        return std::make_shared<Rotation>(
                toDouble(node["rotate"]),
                getOptionalPoint(node, "center", dimensions, {0, 0, 0}),
                toVector3D(node["axis"], Dimensions::Three),
                dimensions);
    }
}

/**
 * Parse a "scale" operator.
 *
 * \param node the node from which to read the operator
 * \param dimensions the expected number of dimensions
 * \return the created operator
 */
OpPtr parseScale(const conduit::Node &node, Dimensions dimensions) {
    verifyObjectFields(node, "scale", FieldSet{}, FieldSet{});
    auto &scaleNode = node["scale"];
    if (scaleNode.dtype().is_number()
        && scaleNode.dtype().number_of_elements() == 1) {
        double factor = scaleNode.to_double();
        return std::make_shared<Scale>(factor, factor, factor, dimensions);
    }
    auto factors = toDoubleVector(scaleNode, toInt(dimensions));
    if (dimensions == Dimensions::Two) {
        factors.emplace_back(1.0);
    }
    return std::make_shared<Scale>(factors[0], factors[1], factors[2],
            dimensions);
}

/**
 * Parse a "matrix" operator.
 *
 * \param node the node from which to read the operator
 * \param dimensions the expected number of dimensions
 * \return the created operator
 */
OpPtr parseMatrix(const conduit::Node &node, Dimensions dimensions) {
    verifyObjectFields(node, "matrix", FieldSet{}, FieldSet{});
    auto &valuesNode = node["matrix"];
    numerics::Matrix<double> matrix = numerics::Matrix<double>::identity(4);
    if (dimensions == Dimensions::Two) {
        auto readEntries = toDoubleVector(valuesNode, 6);
        auto elementIter = readEntries.begin();
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                matrix(i, j) = *elementIter;
                ++elementIter;
            }
            matrix(i, 3) = *elementIter;
            ++elementIter;
        }
    } else {
        auto readEntries = toDoubleVector(valuesNode, 12);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 4; ++j) {
                matrix(i, j) = readEntries[i * 4 + j];
            }
        }
    }
    return std::make_shared<ArbitraryMatrixOperator>(matrix, dimensions);
}

/**
 * Make a slice, ensuring all the values are valid.
 *
 * \param origin the origin of the coordinate system
 * \param normal a vector normal to the plane
 * \param up a vector which defines the positive Y direction
 * \return the created operator
 * \throws std::invalid_argument if any value is invalid
 */
OpPtr makeCheckedSlice(Point3D origin, Vector3D normal, Vector3D up) {
    if (normal.is_zero()) {
        throw std::invalid_argument(
                "The 'normal' vector must not be a zero vector");
    }
    if (!utilities::isNearlyEqual(normal.dot(up), 0.0)) {
        throw std::invalid_argument(
                "The 'normal' and 'up' vectors must be perpendicular");
    }
    return std::make_shared<SliceOperator>(origin, normal, up);
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
        const primal::Vector3D &defaultNormal) {
    double axisIntercept = toDouble(sliceNode[planeName]);

    primal::Point3D defaultOrigin;
    int nonZeroIndex = -1;
    for (int i = 0; i < 3; ++i) {
        defaultOrigin[i] = axisIntercept * defaultNormal[i];
        if (!utilities::isNearlyEqual(defaultNormal[i], 0.0)) {
            nonZeroIndex = i;
        }
    }

    if (!sliceNode.has_child("origin")) {
        return defaultOrigin;
    }

    primal::Point3D givenOrigin = toPoint3D(sliceNode["origin"],
            Dimensions::Three);
    if (givenOrigin[nonZeroIndex] != axisIntercept) {
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
        const primal::Vector3D &defaultNormal) {
    if (!sliceNode.has_child("normal")) {
        return defaultNormal;
    }

    primal::Vector3D givenNormal = toVector3D(sliceNode["normal"],
            Dimensions::Three);
    auto cross = primal::Vector3D::cross_product(
            givenNormal, defaultNormal);
    bool parallel = cross.is_zero();
    if (!parallel) {
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
 * \return the parsed plane
 */
OpPtr readPerpendicularSlice(const conduit::Node &sliceNode,
        char const *planeName, Vector3D const &defaultNormal,
        Vector3D const &defaultUp) {
    verifyObjectFields(sliceNode, planeName, FieldSet{},
            {"origin", "normal", "up"});
    const primal::Vector3D defaultNormalVec{defaultNormal.data()};

    auto origin = getPerpendicularSliceOrigin(sliceNode, planeName,
            defaultNormalVec);
    auto normal = getPerpendicularSliceNormal(sliceNode, defaultNormalVec);
    auto up = getOptionalVector(sliceNode, "up", Dimensions::Three, defaultUp);

    return makeCheckedSlice(origin, normal, up);
}

/**
 * Parse a "slice" operator.
 *
 * \param node the node from which to read the operator
 * \param dimensions the expected number of dimensions
 * \return the created operator
 */
OpPtr parseSlice(const conduit::Node &node, Dimensions dimensions) {
    if (dimensions != Dimensions::Three) {
        throw std::invalid_argument("Cannot do a slice from 2D");
    }
    verifyObjectFields(node, "slice", FieldSet{}, FieldSet{});
    auto &properties = node["slice"];
    if (properties.has_child("x")) {
        return readPerpendicularSlice(properties, "x", {1, 0, 0}, {0, 0, 1});
    } else if (properties.has_child("y")) {
        return readPerpendicularSlice(properties, "y", {0, 1, 0}, {1, 0, 0});
    } else if (properties.has_child("z")) {
        return readPerpendicularSlice(properties, "z", {0, 0, 1}, {0, 1, 0});
    }

    verifyObjectFields(properties, "origin", {"normal", "up"}, FieldSet{});
    return makeCheckedSlice(
            toPoint3D(properties["origin"], Dimensions::Three),
            toVector3D(properties["normal"], Dimensions::Three),
            toVector3D(properties["up"], Dimensions::Three));
}

/**
 * Parse an operator specified via the "ref" command.
 *
 * \param node the node from which to read the operator
 * \param namedOperators a map of named operators from which to get
 * referenced operators
 * \return the created operator
 */
OpPtr parseRef(const conduit::Node &node,
        const NamedOperatorMap &namedOperators) {
    verifyObjectFields(node, "ref", FieldSet{}, FieldSet{});
    std::string const &operatorName = node["ref"].as_string();
    auto opIter = namedOperators.find(operatorName);
    if (opIter == namedOperators.end()) {
        std::string message = "No operator named '";
        message += operatorName;
        message += '\'';
        throw std::invalid_argument(message);
    }
    return opIter->second;
}


/**
 * Parse a single operator
 *
 * \param node the node from which to parse the operator
 * \param dimensions the expected number of dimensions
 * \param namedOperators a map of named operators from which to get
 * referenced operators
 * \return the created operator
 */
OpPtr parseSingleOperator(const conduit::Node &node, Dimensions dimensions,
        const NamedOperatorMap &namedOperators) {

    std::unordered_map<std::string, OperatorParser> parsers{
            {"translate", parseTranslate},
            {"rotate", parseRotate},
            {"scale", parseScale},
            {"matrix", parseMatrix},
            {"slice", parseSlice},
            {"ref", [&namedOperators](const conduit::Node &node, Dimensions) {
                return parseRef(node, namedOperators);
            }}
    };

    for (auto &entry : parsers) {
        if (node.has_child(entry.first)) {
            return entry.second(node, dimensions);
        }
    }

    std::string message = "Invalid transformation: \n";
    message += node.to_json_default();
    throw std::invalid_argument(message);
}

}

std::shared_ptr<const GeometryOperator> parseGeometryOperators(
        const conduit::Node &node, Dimensions initialDimensions,
        const NamedOperatorMap &namedOperators) {
    auto composite = std::make_shared<CompositeOperator>();
    Dimensions currentDimensions = initialDimensions;
    for (auto iter = node.children(); iter.has_next();) {
        auto op = parseSingleOperator(iter.next(), currentDimensions,
                namedOperators);
        composite->addOperator(op);
        currentDimensions = composite->endDims();
    }
    return composite;
}

NamedOperatorMap parseNamedGeometryOperators(const conduit::Node &node,
        Dimensions initialDimensions) {
    NamedOperatorMap namedOperators;

    for (auto iter = node.children(); iter.has_next();) {
        auto &namedOperator = iter.next();
        std::string name = namedOperator["name"].as_string();

        Dimensions dimensions = initialDimensions;
        if (namedOperator.has_child("initial_dimensions")) {
            dimensions = toDimensions(namedOperator["initial_dimensions"]);
        }

        auto op = parseGeometryOperators(namedOperator["value"],
                dimensions, namedOperators);
        namedOperators.insert({name, op});
    }
    return namedOperators;
}

}}}
