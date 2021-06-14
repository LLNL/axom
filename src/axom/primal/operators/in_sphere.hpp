// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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
#include "axom/core/numerics/Determinants.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Tests whether a query point lies inside a 2D triangle's circumcircle
 *
 * A triangle's circumscircle is the unique circle (i.e. a 2-sphere) that
 * passes through each of its three vertices.
 *
 * \param [in] q the query point
 * \param [in] p0 the first vertex of the triangle
 * \param [in] p1 the second vertex of the triangle
 * \param [in] p2 the third vertex of the triangle
 * \return true if the point is inside the circumcircle, false if it is on
 * the circle's boundary or outside the circle
 */
template <typename T>
inline bool in_sphere(const Point<T, 2>& q,
                      const Point<T, 2>& p0,
                      const Point<T, 2>& p1,
                      const Point<T, 2>& p2)
{
  // clang-format off
  double det = axom::numerics::determinant(
    1.0, p0[0], p0[1], p0[0]*p0[0] + p0[1]*p0[1],
    1.0, p1[0], p1[1], p1[0]*p1[0] + p1[1]*p1[1],
    1.0, p2[0], p2[1], p2[0]*p2[0] + p2[1]*p2[1],
    1.0,  q[0],  q[1],  q[0]* q[0] +  q[1]* q[1]);
  // clang-format on

  return det < 0;
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
 * \return true if the point is inside the circumsphere, false if it is on
 * the sphere's boundary or outside the sphere
 */
template <typename T>
inline bool in_sphere(const Point<T, 3>& q,
                      const Point<T, 3>& p0,
                      const Point<T, 3>& p1,
                      const Point<T, 3>& p2,
                      const Point<T, 3>& p3)
{
  double mat_val[] = {
    1.0, p0[0], p0[1], p0[2], p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2],
    1.0, p1[0], p1[1], p1[2], p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2],
    1.0, p2[0], p2[1], p2[2], p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2],
    1.0, p3[0], p3[1], p3[2], p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2],
    1.0, q[0],  q[1],  q[2],  q[0] * q[0] + q[1] * q[1] + q[2] * q[2]};

  axom::numerics::Matrix<double> mat(5, 5, mat_val, true);

  double det = axom::numerics::determinant(mat);
  return det < 0;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_IN_SPHERE_H_
