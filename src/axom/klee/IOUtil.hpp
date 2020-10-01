// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_IO_UTIL_HPP
#define AXOM_KLEE_IO_UTIL_HPP

#include <vector>

#include "conduit.hpp"

#include "axom/klee/Dimensions.hpp"

namespace axom { namespace klee { namespace internal {

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

}}}


#endif //AXOM_KLEE_IO_UTIL_HPP
