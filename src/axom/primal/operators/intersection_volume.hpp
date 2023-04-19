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
 * \param [in] checkSign If true (default is false), checks the volumes of the
 *             shapes are positive. If volume is negative, order of some
 *             vertices will be swapped.
 *
 * \return Intersection volume between the hexahedron and tetrahedron
 *
 * \note checkSign flag does not guarantee the shapes' vertex orders
 *       will be valid. It is the responsiblity of the caller to pass
 *       shapes with a valid vertex order.
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Hexahedron<T, 3>& hex,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10,
                                       bool checkSign = false)
{
  return clip(hex, tet, eps, checkSign).volume();
}

/*!
 * \brief Finds the intersection volume between a tetrahedron and a hexahedron
 *
 * \param [in] hex The tetrahedron
 * \param [in] tet The hexahedron
 * \param [in] eps The tolerance for determining the intersection
 * \param [in] checkSign If true (default is false), checks the volumes of the
 *             shapes are positive. If volume is negative, order of some
 *             vertices will be swapped.
 *
 * \return Intersection volume between the tetrahedron and hexahedron
 *
 * \note checkSign flag does not guarantee the shapes' vertex orders
 *       will be valid. It is the responsiblity of the caller to pass
 *       shapes with a valid vertex order.
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Tetrahedron<T, 3>& tet,
                                       const Hexahedron<T, 3>& hex,
                                       double eps = 1.e-10,
                                       bool checkSign = false)
{
  return intersection_volume(hex, tet, eps, checkSign);
}

/*!
 * \brief Finds the intersection volume between a octahedron and a tetrahedron
 *
 * \param [in] oct The octahedron
 * \param [in] tet The tetrahedron
 * \param [in] eps The tolerance for determining the intersection
 * \param [in] checkSign If true (default is false), checks the volumes of the
 *             shapes are positive. If volume is negative, order of some
 *             vertices will be swapped.
 *
 * \return Intersection volume between the octahedron and tetrahedron
 *
 * \note checkSign flag does not guarantee the shapes' vertex orders
 *       will be valid. It is the responsiblity of the caller to pass
 *       shapes with a valid vertex order.
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Octahedron<T, 3>& oct,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10,
                                       bool checkSign = false)
{
  return clip(oct, tet, eps, checkSign).volume();
}

/*!
 * \brief Finds the intersection volume between a tetrahedron and a octahedron
 *
 * \param [in] oct The tetrahedron
 * \param [in] tet The octahedron
 * \param [in] eps The tolerance for determining the intersection
 * \param [in] checkSign If true (default is false), checks the volumes of the
 *             shapes are positive. If volume is negative, order of some
 *             vertices will be swapped.
 *
 * \return Intersection volume between the tetrahedron and octahedron
 *
 * \note checkSign flag does not guarantee the shapes' vertex orders
 *       will be valid. It is the responsiblity of the caller to pass
 *       shapes with a valid vertex order.
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Tetrahedron<T, 3>& tet,
                                       const Octahedron<T, 3>& oct,
                                       double eps = 1.e-10,
                                       bool checkSign = false)
{
  return intersection_volume(oct, tet, eps, checkSign);
}

/*!
 * \brief Finds the intersection volume between a tetrahedron and another
 *        tetrahedron
 *
 * \param [in] tet1 The tetrahedron
 * \param [in] tet2 The other tetrahedron
 * \param [in] eps The tolerance for determining the intersection
 * \param [in] checkSign If true (default is false), checks the volumes of the
 *             shapes are positive. If volume is negative, order of some
 *             vertices will be swapped.
 *
 * \return Intersection volume between the tetrahedra
 *
 * \note checkSign flag does not guarantee the shapes' vertex orders
 *       will be valid. It is the responsiblity of the caller to pass
 *       shapes with a valid vertex order.
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Tetrahedron<T, 3>& tet1,
                                       const Tetrahedron<T, 3>& tet2,
                                       double eps = 1.e-10,
                                       bool checkSign = false)
{
  return clip(tet1, tet2, eps, checkSign).volume();
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECTION_VOLUME_HPP_
