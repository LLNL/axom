// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersection_volume.hpp
 *
 * \brief Consists of a set of methods to find the volume of
 *        intersection (clipping) between a primal primitive and
 *        another primal primitive
 */

#ifndef AXOM_PRIMAL_INTERSECTION_VOLUME_HPP_
#define AXOM_PRIMAL_INTERSECTION_VOLUME_HPP_

#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "axom/primal/operators/clip.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Finds the intersection volume between a hexahedron and a tetrahedron
 *
 * \param [in] hex The hexahedron
 * \param [in] tet The tetrahedron
 * \param [in] eps The tolerance for determining the intersection
 * \return Intersection volume between the hexahedron and tetrahedron
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Hexahedron<T, 3>& hex,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10)
{
  return clip(hex, tet, eps).volume();
}

/*!
 * \brief Finds the intersection volume between a octahedron and a tetrahedron
 *
 * \param [in] oct The octahedron
 * \param [in] tet The tetrahedron
 * \param [in] eps The tolerance for determining the intersection
 * \return Intersection volume between the octahedron and tetrahedron
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Octahedron<T, 3>& oct,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10)
{
  return clip(oct, tet, eps).volume();
}

/*!
 * \brief Finds the intersection volume between a tetrahedron and another
 *        tetrahedron
 *
 * \param [in] tet1 The tetrahedron
 * \param [in] tet2 The other tetrahedron
 * \param [in] eps The tolerance for determining the intersection
 * \return Intersection volume between the tetrahedra
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Tetrahedron<T, 3>& tet1,
                                       const Tetrahedron<T, 3>& tet2,
                                       double eps = 1.e-10)
{
  return clip(tet1, tet2, eps).volume();
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECTION_VOLUME_HPP_
