// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_GEOMETRYOPERATORSIO_HPP
#define AXOM_KLEE_GEOMETRYOPERATORSIO_HPP

#include <memory>

#include "conduit.hpp"

#include "axom/klee/Dimensions.hpp"

namespace axom { namespace klee {

class GeometryOperator;

namespace internal {

/**
 * Parse a GeometryOperator from the given node.
 *
 * \param node a conduit node representing the operator. This should be a list
 * where each entry is its own operator.
 * \param initialDimensions the number of dimensions in which the first
 * operator must operate
 * \return A GeometryOperator describing the parsed values.
 */
std::shared_ptr<const GeometryOperator> parseGeometryOperators(
        const conduit::Node &node, Dimensions initialDimensions);

}}}

#endif //AXOM_KLEE_GEOMETRYOPERATORSIO_HPP
