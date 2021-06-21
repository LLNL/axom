// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file closest_point.hpp
 *
 * \brief Consists of a set of methods that compute the closest point on a
 *  geometric primitive B from another geometric primitive A.
 *
 */

#ifndef AXOM_PRIMAL_CLOSEST_POINT_HPP_
#define AXOM_PRIMAL_CLOSEST_POINT_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Computes the closest point from a point, P, to a given triangle.
 *
 * \param [in] P the query point
 * \param [in] tri user-supplied triangle.
 * \param [out] loc int pointer to store location of closest point (optional).
 * \return cp the closest point from a point P and a triangle.
 *
 * \note If the optional int pointer is supplied for `loc`, the method returns
 *  the location of the closest point, which is illustrated in the schematic
 *  diagram below and encoded as follows:
 * <ul>
 *  <li> loc \f$ \in [0,2] \f$, loc corresponds to the triangle node index </li>
 *  <li> loc \f$ \in [-3,-1] \f$, abs(loc) corresponds to an edge </li>
 *  <li> loc >= 3, loc is on a triangle face </li>
 * </ul>
 *
 * \verbatim
 *
 *            2
 *           /\
 *    (-3)--/  \--(-2)
 *         /    \
 *        /_ _ _ \
 *       0   |    1
 *           |
 *         (-1)
 *
 * \endverbatim
 *
 * \pre NDIMS==2 || NDIMS==3
 *
 * \note Implementation is based on "Real Time Collision Detection,
 *  Chapter 5.1.5 Closest Point on Triangle to Point".
 */
template <typename T, int NDIMS>
inline Point<T, NDIMS> closest_point(const Point<T, NDIMS>& P,
                                     const Triangle<T, NDIMS>& tri,
                                     int* loc = nullptr)
{
// convenience macros to access triangle vertices
#define A(t) t[0]
#define B(t) t[1]
#define C(t) t[2]

  // Check if P in vertex region outside A
  Vector<T, NDIMS> ab(A(tri), B(tri));
  Vector<T, NDIMS> ac(A(tri), C(tri));
  Vector<T, NDIMS> ap(A(tri), P);
  T d1 = Vector<T, NDIMS>::dot_product(ab, ap);
  T d2 = Vector<T, NDIMS>::dot_product(ac, ap);
  if(d1 <= 0.0f && d2 <= 0.0f)
  {
    // A is the closest point
    if(loc != nullptr)
    {
      *loc = 0;
    }

    return (A(tri));

  }  // END if

  //----------------------------------------------------------------------------
  // Check if P in vertex region outside B
  Vector<T, NDIMS> bp(B(tri), P);
  T d3 = Vector<T, NDIMS>::dot_product(ab, bp);
  T d4 = Vector<T, NDIMS>::dot_product(ac, bp);
  if(d3 >= 0.0f && d4 <= d3)
  {
    // B is the closest point
    if(loc != nullptr)
    {
      *loc = 1;
    }

    return (B(tri));

  }  // END if

  //----------------------------------------------------------------------------
  // Check if P in edge region of AB
  T vc = d1 * d4 - d3 * d2;
  if(vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
  {
    T v = d1 / (d1 - d3);
    Vector<T, NDIMS> v_ab = ab * v;

    double x = A(tri)[0] + v_ab[0];
    double y = A(tri)[1] + v_ab[1];
    double z = (NDIMS == 3) ? A(tri)[2] + v_ab[2] : 0.0;

    if(loc != nullptr)
    {
      *loc = -1;
    }

    return (Point<T, NDIMS>::make_point(x, y, z));
  }  // END if

  //----------------------------------------------------------------------------
  // Check if P in vertex region outside C
  Vector<T, NDIMS> cp(C(tri), P);
  T d5 = Vector<T, NDIMS>::dot_product(ab, cp);
  T d6 = Vector<T, NDIMS>::dot_product(ac, cp);
  if(d6 >= 0.0f && d5 <= d6)
  {
    // C is the closest point
    if(loc != nullptr)
    {
      *loc = 2;
    }

    return (C(tri));
  }

  //----------------------------------------------------------------------------
  // Check if P in edge region of AC
  T vb = d5 * d2 - d1 * d6;
  if(vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
  {
    T w = d2 / (d2 - d6);
    Vector<T, NDIMS> w_ac = ac * w;

    double x = A(tri)[0] + w_ac[0];
    double y = A(tri)[1] + w_ac[1];
    double z = (NDIMS == 3) ? A(tri)[2] + w_ac[2] : 0.0;

    if(loc != nullptr)
    {
      *loc = -3;
    }

    return (Point<T, NDIMS>::make_point(x, y, z));
  }  // END if

  //----------------------------------------------------------------------------
  // Check if P in edge region of BC
  T va = d3 * d6 - d5 * d4;
  if(va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)
  {
    T w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    Vector<T, NDIMS> bc(B(tri), C(tri));
    Vector<T, NDIMS> w_bc = bc * w;

    double x = B(tri)[0] + w_bc[0];
    double y = B(tri)[1] + w_bc[1];
    double z = (NDIMS == 3) ? B(tri)[2] + w_bc[2] : 0.0;

    if(loc != nullptr)
    {
      *loc = -2;
    }

    return (Point<T, NDIMS>::make_point(x, y, z));
  }  // END if

  //----------------------------------------------------------------------------
  // P is inside face region
  T denom = 1.0f / (va + vb + vc);
  T v = vb * denom;
  T w = vc * denom;
  Vector<T, NDIMS> N = (ab * v) + (ac * w);

  double x = A(tri)[0] + N[0];
  double y = A(tri)[1] + N[1];
  double z = (NDIMS == 3) ? A(tri)[2] + N[2] : 0.0;

  if(loc != nullptr)
  {
    *loc = Triangle<T, NDIMS>::NUM_TRI_VERTS;
  }

  return (Point<T, NDIMS>::make_point(x, y, z));

#undef A
#undef B
#undef C
}

/*!
 * \brief Computes the closest point from a point to a given OBB.
 *
 * \param [in] pt the query pt.
 * \param [in] obb user-supplied oriented bounding box.
 * \return cp the closest point from a point pt and an OBB.
 */
template <typename T, int NDIMS>
inline Point<T, NDIMS> closest_point(const Point<T, NDIMS>& pt,
                                     const OrientedBoundingBox<T, NDIMS>& obb)
{
  Vector<T, NDIMS> e = obb.getExtents();
  const Vector<T, NDIMS>* u = obb.getAxes();

  Vector<T, NDIMS> pt_l = obb.toLocal(pt);
  Vector<T, NDIMS> res(obb.getCentroid());

  for(int i = 0; i < NDIMS; i++)
  {
    // since the local coordinates are individually constrained, we can simply
    // choose the "best" local coordinate in each axis direction
    if(pt_l[i] <= e[i] && pt_l[i] >= -e[i])
    {
      res += pt_l[i] * u[i];
    }
    else if(pt_l[i] > e[i])
    {
      res += e[i] * u[i];
    }
    else
    {
      res -= e[i] * u[i];
    }
  }

  return Point<T, NDIMS>(res.array());
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CLOSEST_POINT_HPP_
