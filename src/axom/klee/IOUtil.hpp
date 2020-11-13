// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_IO_UTIL_HPP
#define AXOM_KLEE_IO_UTIL_HPP

#include <tuple>
#include <vector>

#include "conduit.hpp"

#include "axom/klee/Dimensions.hpp"
#include "axom/klee/Units.hpp"

namespace axom
{
namespace klee
{
namespace internal
{
/**
 * Convert a Conduit node to a vector of doubles.
 *
 * \param listNode the Conduit node containing the list of values
 * \param expectedSize the number of entries that should be in the list
 * \return a vector with the values in the list
 * \throws std::invalid_argument if the node does not represent a list
 * of exactly the specified number of elements.
 */
std::vector<double> toDoubleVector(const conduit::Node &listNode,
                                   std::size_t expectedSize);

/**
 * Convert a Conduit node to a double, or throw an exception with a nice
 * error message when this is not possible.
 *
 * \param value the value as a Conduit node.
 * \return the value of the node as a double
 * \throws std::invalid_argument if the node is not a double
 */
double toDouble(const conduit::Node &value);

/**
 * Convert a Conduit node to Dimensions. Must be an integer with a value of
 * 2 or 3/
 * \param dimensionsNode the node to convert
 * \return the number of dimensions
 * \throws std::invalid_argument if the node is not a valid dimension value
 */
Dimensions toDimensions(const conduit::Node &dimensionsNode);

/**
 * Get the start and end units in a node.
 *
 * The node may either have a "units" field, or a "start_units" and "end_units".
 * In the first case, "units" will be used for both the start and end. In the
 * second, both must be present. In the case where no units are present at
 * all, both returned units will be LengthUnit::unspecified.
 *
 * \param node the node from which to get the units
 * \return the start and end units
 * \throws std::invalid_argument if an invalid combination of fields is specified
 */
std::tuple<LengthUnit, LengthUnit> getOptionalStartAndEndUnits(
  const conduit::Node &node);

/**
 * Get the start and end units in a node.
 *
 * The node may either have a "units" field, or a "start_units" and "end_units".
 * In the first case, "units" will be used for both the start and end. In the
 * second, both must be present.
 *
 * \param node the node from which to get the units
 * \return the start and end units
 * \throws std::invalid_argument if an invalid combination of fields is
 * specified or if no units are specified.
 */
std::tuple<LengthUnit, LengthUnit> getStartAndEndUnits(const conduit::Node &node);

}  // namespace internal
}  // namespace klee
}  // namespace axom

#endif  //AXOM_KLEE_IO_UTIL_HPP
