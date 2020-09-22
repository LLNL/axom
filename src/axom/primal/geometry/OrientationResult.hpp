// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_ORIENTATIONRESULT_HPP_
#define PRIMAL_ORIENTATIONRESULT_HPP_

/*!
 * \file
 *
 * \brief Defines the OrientationResult enum which defines possible return
 *  values for orientation tests in between different geometric primitives
 */

namespace axom
{
namespace primal
{
/*!
 * \brief Enumerates possible return values for orientation tests.
 */
enum OrientationResult
{
  ON_BOUNDARY,      /*!< primitive is on the boundary of a primitive      */
  ON_POSITIVE_SIDE, /*!< primitive is on the positive side of a primitive */
  ON_NEGATIVE_SIDE  /*!< primitive is on the negative side of a primitive */
};

}  // namespace primal
}  // namespace axom

#endif /* PRIMAL_ORIENTATIONRESULT_HPP_ */
