// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file winding_number.hpp
 *
 * \brief Consists of methods to compute winding numbers for points 
 *        with respect to various geometric objects.
 */

#ifndef AXOM_PRIMAL_WINDING_NUMBER_HPP_
#define AXOM_PRIMAL_WINDING_NUMBER_HPP_

// Axom includes
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/detail/winding_number_impl.hpp"

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
/*
 * \brief Compute the winding number with respect to a line segment
 *
 * \param [in] q The query point to test
 * \param [in] c0 The initial point of the line segment
 * \param [in] c1 The terminal point of the line segment
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * The winding number for a point with respect to a straight line
 * is the signed angle subtended by the query point to each endpoint.
 *
 * \return double The winding number
 */
template <typename T>
double winding_number(const Point<T, 2>& q, const Segment<T, 2>& s, double edge_tol)
{
  return linear_winding_number(q, s[0], s[1], edge_tol);
}

/*!
 * \brief Computes the winding number for a point and a polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [in] useStrictInclusion If true, points on the boundary are considered exterior.
 * \param [in] EPS The tolerance level for collinearity
 * 
 * Uses an adapted ray-casting approach that counts quarter-rotation
 * of vertices around the query point. 
 *
 * Directly uses algorithm in 
 * Kai Hormann, Alexander Agathos, "The point in polygon problem for arbitrary polygons"
 * Computational Geometry, Volume 20, Issue 3, 2001,
 * 
 * \return The integer winding number
 */
template <typename T>
int winding_number(const Point<T, 2>& R,
                   const Polygon<T, 2>& P,
                   bool useStrictInclusion = false,
                   double EPS = 1e-8)
{
  const int nverts = P.numVertices();

  // If the query is a vertex, return a value interpreted
  //  as "inside" by evenodd or nonzero protocols
  if(axom::utilities::isNearlyEqual(P[0][0], R[0], EPS) &&
     axom::utilities::isNearlyEqual(P[0][1], R[1], EPS))
    return !useStrictInclusion;

  int winding_num = 0;
  for(int i = 0; i < nverts; i++)
  {
    int j = (i == nverts - 1) ? 0 : i + 1;

    if(axom::utilities::isNearlyEqual(P[j][1], R[1], EPS))
    {
      if(axom::utilities::isNearlyEqual(P[j][0], R[0], EPS))
        return !useStrictInclusion;  // On vertex
      else if(P[i][1] == R[1] && ((P[j][0] > R[0]) == (P[i][0] < R[0])))
        return !useStrictInclusion;  // On horizontal edge
    }

    // Check if edge crosses horizontal line
    if((P[i][1] < R[1]) != (P[j][1] < R[1]))
    {
      double det;
      if(P[i][0] >= R[0])
      {
        if(P[j][0] > R[0])
          winding_num += 2 * (P[j][1] > P[i][1]) - 1;
        else
        {
          // clang-format off
          det = axom::numerics::determinant(P[i][0] - R[0], P[j][0] - R[0],
                                            P[i][1] - R[1], P[j][1] - R[1]);
          // clang-format on

          // On edge
          if(axom::utilities::isNearlyEqual(det, 0.0, EPS))
            return !useStrictInclusion;

          // Check if edge intersects horitonal ray to the right of R
          if((det > 0) == (P[j][1] > P[i][1]))
            winding_num += 2 * (P[j][1] > P[i][1]) - 1;
        }
      }
      else
      {
        if(P[j][0] > R[0])
        {
          // clang-format off
          det = axom::numerics::determinant(P[i][0] - R[0], P[j][0] - R[0],
                                            P[i][1] - R[1], P[j][1] - R[1]);
          // clang-format on

          // On edge
          if(axom::utilities::isNearlyEqual(det, 0.0, EPS))
            return !useStrictInclusion;

          // Check if edge intersects horitonal ray to the right of R
          if((det > 0) == (P[j][1] > P[i][1]))
            winding_num += 2 * (P[j][1] > P[i][1]) - 1;
        }
      }
    }
  }

  return winding_num;
}

/*!
 * \brief Computes the generalized winding number for a single Bezier curve
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The Bezier curve object 
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the winding number using a recursive, bisection algorithm,
 * using nearly-linear Bezier curves as a base case.
 * 
 * \return float the generalized winding number.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const BezierCurve<T, 2>& c,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  return detail::curve_winding_number_recursive(q, c, false, edge_tol, EPS);
}

/*!
 * \brief Computes the generalized winding number for a curved polygon
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The CurvedPolygon object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the winding number using a recursive, bisection algorithm.
 * Iterates over the edges of a curved polygon object, and uses nearly-linear 
 * Bezier curves as a base case.
 * 
 * \return float the generalized winding number.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const CurvedPolygon<T, 2>& cpoly,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  double ret_val = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
    ret_val +=
      detail::curve_winding_number_recursive(q, cpoly[i], false, edge_tol, EPS);

  return ret_val;
}

/*!
 * \brief Computes the solid angle winding number for a 3D triangle
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The Triangle object
 *
 * Computes the winding number using the formula from [Barill 2018]
 * 
 * \return float the generalized winding number.
 */
template <typename T>
double winding_number(const Point<T, 3>& q, const Triangle<T, 3>& tri)
{
  if(tri.area() == 0) return 0;

  Vector<T, 3> a(q, tri[0]), b(q, tri[1]), c(q, tri[2]);

  // Compute norms. Possibly return early
  double a_norm = a.norm(), b_norm = b.norm(), c_norm = c.norm();
  if(a_norm == 0 || b_norm == 0 || c_norm == 0) return 0;

  double num = axom::utilities::abs(Vector<T, 3>::scalar_triple_product(a, b, c));
  double denom = a_norm * b_norm * c_norm +
    a_norm * Vector<T, 3>::dot_product(b, c) +
    b_norm * Vector<T, 3>::dot_product(a, c) +
    c_norm * Vector<T, 3>::dot_product(a, b);

  return 2 * atan(num / denom);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_H_
