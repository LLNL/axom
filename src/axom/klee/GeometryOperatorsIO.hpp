// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_GEOMETRYOPERATORSIO_HPP
#define AXOM_KLEE_GEOMETRYOPERATORSIO_HPP

#include <memory>
#include <string>
#include <unordered_map>

#include "conduit.hpp"

#include "axom/klee/Dimensions.hpp"

namespace axom
{
namespace klee
{
class GeometryOperator;

namespace internal
{
using NamedOperatorMap =
  std::unordered_map<std::string, std::shared_ptr<const GeometryOperator>>;

/**
 * Parse a GeometryOperator from the given node.
 *
 * \param node a conduit node representing the operator. This should be a list
 * where each entry is its own operator.
 * \param initialDimensions the number of dimensions in which the first
 * operator must operate
 * \param namedOperators a map of named operators which contains operators
 * which may be referenced via the "ref" syntax.
 * \return A GeometryOperator describing the parsed values.
 */
std::shared_ptr<const GeometryOperator> parseGeometryOperators(
  const conduit::Node &node,
  Dimensions initialDimensions,
  const NamedOperatorMap &namedOperators);

/**
 * Parse named geometry operators.
 *
 * \param node the node from which to parse the operators. This should be
 * the value of the "named_operators" node in the specification.
 * \param initialDimensions the number of dimensions that operators should
 * start in, unless they specify otherwise.
 * \return a map of names to operators
 */
NamedOperatorMap parseNamedGeometryOperators(const conduit::Node &node,
                                             Dimensions initialDimensions);

}  // namespace internal
}  // namespace klee
}  // namespace axom

#endif  //AXOM_KLEE_GEOMETRYOPERATORSIO_HPP
