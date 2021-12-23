// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file in_sphere.hpp
 *
 * \brief Consists of methods that test whether a given query point is
 * inside the unique sphere circumscribing a 2D triangle or a 3D tetrahehdron.
 *
 * This is a well known computational geometry primitive.  For reference,
 * see Section 3.1.6.4 in "Real-time collision detection" by C. Ericson.
 */

#ifndef AXOM_PRIMAL_IN_SPHERE_H_
#define AXOM_PRIMAL_IN_SPHERE_H_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/OrientationResult.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Tests whether a query point lies inside a 2D triangle's circumcircle
 *
 * A triangle's circumcircle is the unique circle (i.e. a 2-sphere) that
 * passes through each of its three vertices.
 *
 * \param [in] q the query point
 * \param [in] p0 the first vertex of the triangle
 * \param [in] p1 the second vertex of the triangle
 * \param [in] p2 the third vertex of the triangle
 * \param [in] EPS tolerance for determining if \a q is on the boundary. Default: 1e-8.
 * \return true if the point is inside the circumcircle, false if it is on
 * the circle's boundary or outside the circle
 */
template <typename T>
inline bool in_sphere(const Point<T, 2>& q,
                      const Point<T, 2>& p0,
                      const Point<T, 2>& p1,
                      const Point<T, 2>& p2,
                      double EPS = 1e-8)
{
  return in_sphere(q, Triangle<T, 2>(p0, p1, p2), EPS);
}

/*!
 * \brief Tests whether a query point lies inside a 2D triangle's circumcircle
 *
 * \param [in] q the query point
 * \param [in] tri the triangle
 * \param [in] EPS tolerance for determining if \a q is on the boundary. Default: 1e-8.
 * \see in_sphere
 */
template <typename T>
inline bool in_sphere(const Point<T, 2>& q,
                      const Triangle<T, 2>& tri,
                      double EPS = 1e-8)
{
  const auto sphere = tri.circumsphere();
  return sphere.getOrientation(q.data(), EPS) == ON_NEGATIVE_SIDE;
}

/*!
 * \brief Tests whether a query point lies inside a 3D tetrahedron's
 * circumsphere
 *
 * A tetrahedron's circumsphere is the unique sphere that passes through each
 * of its four vertices.
 *
 * \param [in] q the query point
 * \param [in] p0 the first vertex of the tetrahedron
 * \param [in] p1 the second vertex of the tetrahedron
 * \param [in] p2 the third vertex of the tetrahedron
 * \param [in] p3 the fourth vertex of the tetrahedron
 * \param [in] EPS tolerance for determining if \a q is on the boundary. Default: 1e-8.
 * \return true if the point is inside the circumsphere, false if it is on
 * the sphere's boundary or outside the sphere
 */
template <typename T>
inline bool in_sphere(const Point<T, 3>& q,
                      const Point<T, 3>& p0,
                      const Point<T, 3>& p1,
                      const Point<T, 3>& p2,
                      const Point<T, 3>& p3,
                      double EPS = 1e-8)
{
  return in_sphere(q, Tetrahedron<T, 3>(p0, p1, p2, p3), EPS);
}

/*!
 * \brief Tests whether a query point lies inside a 3D tetrahedron's circumsphere
 *
 * \param [in] q the query point
 * \param [in] tet the tetrahedron
 * \param [in] EPS tolerance for determining if \a q is on the boundary. Default: 1e-8.
 * \see in_sphere
 */
template <typename T>
inline bool in_sphere(const Point<T, 3>& q,
                      const Tetrahedron<T, 3>& tet,
                      double EPS = 1e-8)
{
  const auto sphere = tet.circumsphere();
  return sphere.getOrientation(q.data(), EPS) == ON_NEGATIVE_SIDE;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_IN_SPHERE_H_
