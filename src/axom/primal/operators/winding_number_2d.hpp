// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file winding_number_2d.hpp
 *
 * \brief Consists of methods to compute the generalized winding number (GWN) 
 *        for points with respect to various 2D geometric objects.
 */

#ifndef AXOM_PRIMAL_WINDING_NUMBER_2D_HPP_
#define AXOM_PRIMAL_WINDING_NUMBER_2D_HPP_

// Axom includes
#include "axom/core.hpp"
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/detail/winding_number_2d_impl.hpp"

// C++ includes
#include <cmath>

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

namespace axom
{
namespace primal
{
/*
 * \brief Compute the GWN for a 2D point wrt a 2D line segment
 *
 * \param [in] q The query point to test
 * \param [in] s The line segment
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * \return The GWN
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const Segment<T, 2>& s,
                      double edge_tol = 1e-8)
{
  return detail::linear_winding_number(q, s[0], s[1], edge_tol);
}

/*
 * \brief Compute the winding number for a 2D point wrt a 2D triangle
 *
 * \param [in] q The query point to test
 * \param [in] tri The triangle
 * \param [in] includeBoundary If true, points on the boundary are considered interior.
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * The triangle is assumed to be closed, so the winding number is an integer
 * 
 * \return The integer winding number
 */
template <typename T>
int winding_number(const Point<T, 2>& q,
                   const Triangle<T, 2>& tri,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8)
{
  return winding_number(
    q,
    Polygon<T, 2>(axom::Array<Point<T, 2>>({tri[0], tri[1], tri[2]})),
    includeBoundary,
    edge_tol);
}

/*!
 * \brief Computes the winding number for a 2D point wrt a 2D polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [in] includeBoundary If true, points on the boundary are considered interior
 * \param [in] isOnEdge An optional return parameter if the point is on the boundary
 * \param [in] edge_tol The distance at which a point is considered on the boundary
 * 
 * Uses an adapted ray-casting approach that counts quarter-rotation
 * of vertices around the query point. Current policy is to return 1 on edges
 * without strict inclusion, 0 on edges with strict inclusion.
 *
 * The polygon is assumed to be closed, so the winding number is an integer
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
                   bool& isOnEdge,
                   bool includeBoundary,
                   double edge_tol)
{
  const int nverts = P.numVertices();
  const double edge_tol_2 = edge_tol * edge_tol;
  isOnEdge = false;

  int winding_num = 0;
  for(int i = 0; i < nverts; i++)
  {
    int j = (i == nverts - 1) ? 0 : i + 1;

    // Check if the point is on the edge up to some tolerance
    if(squared_distance(R, Segment<T, 2>(P[i], P[j])) <= edge_tol_2)
    {
      isOnEdge = true;
      return includeBoundary ? 1 : 0;
    }

    // Check if edge crosses horizontal line
    if((P[i][1] < R[1]) != (P[j][1] < R[1]))
    {
      if(P[i][0] >= R[0])
      {
        if(P[j][0] > R[0])
        {
          winding_num += 2 * (P[j][1] > P[i][1]) - 1;
        }
        else
        {
          // clang-format off
          double det = axom::numerics::determinant(P[i][0] - R[0], P[j][0] - R[0],
                                                   P[i][1] - R[1], P[j][1] - R[1]);
          // clang-format on

          // Check if edge intersects horitonal ray to the right of R
          if((det > 0) == (P[j][1] > P[i][1]))
          {
            winding_num += 2 * (P[j][1] > P[i][1]) - 1;
          }
        }
      }
      else
      {
        if(P[j][0] > R[0])
        {
          // clang-format off
          double det = axom::numerics::determinant(P[i][0] - R[0], P[j][0] - R[0],
                                                   P[i][1] - R[1], P[j][1] - R[1]);
          // clang-format on

          // Check if edge intersects horitonal ray to the right of R
          if((det > 0) == (P[j][1] > P[i][1]))
          {
            winding_num += 2 * (P[j][1] > P[i][1]) - 1;
          }
        }
      }
    }
  }

  return winding_num;
}

/*!
 * \brief Computes the winding number for a 2D point wrt a 2D polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [in] includeBoundary If true, points on the boundary are considered interior
 * \param [in] edge_tol The distance at which a point is considered on the boundary
 * 
 * Computes the integer winding number for a polygon without an additional
 *  return parameter for whether the point is on the boundary.
 * 
 * \return The integer winding number
 */
template <typename T>
int winding_number(const Point<T, 2>& R,
                   const Polygon<T, 2>& P,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8)
{
  bool isOnEdge = false;
  return winding_number(R, P, isOnEdge, includeBoundary, edge_tol);
}

/*!
 * \brief Computes the GWN for a 2D point wrt a 2D Bezier curve
 *
 * \param [in] query The query point to test
 * \param [in] c The Bezier curve object 
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN using a recursive, bisection algorithm
 * that constructs a polygon with the same *integer* WN as 
 * the curve closed with a linear segment. The *generalized* WN 
 * of the closing line is then subtracted from the integer WN to
 * return the GWN of the original curve.
 *  
 * Nearly-linear Bezier curves are the base case for recursion.
 * 
 * See Algorithm 2 in
 *  Jacob Spainhour, David Gunderman, and Kenneth Weiss. 2024. 
 *  Robust Containment Queries over Collections of Rational Parametric Curves via Generalized Winding Numbers. 
 *  ACM Trans. Graph. 43, 4, Article 38 (July 2024)
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const BezierCurve<T, 2>& c,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  const int ord = c.getOrder();
  if(ord <= 0) return 0.0;

  // Early return is possible for most points + curves
  if(!c.boundingBox().expand(edge_tol).contains(q))
  {
    return detail::linear_winding_number(q, c[0], c[ord], edge_tol);
  }

  // The first vertex of the polygon is the t=0 point of the curve
  Polygon<T, 2> approximating_polygon(1);
  approximating_polygon.addVertex(c[0]);

  // Need to keep a running total of the GWN to account for
  //  the winding number of coincident points
  double gwn = 0.0;
  bool isCoincident = false;
  detail::construct_approximating_polygon(q,
                                          c,
                                          false,
                                          edge_tol,
                                          EPS,
                                          approximating_polygon,
                                          gwn,
                                          isCoincident);

  // The last vertex of the polygon is the t=1 point of the curve
  approximating_polygon.addVertex(c[ord]);

  // Compute the integer winding number of the closed curve
  bool isOnEdge = false;
  double closed_curve_wn =
    winding_number(q, approximating_polygon, isOnEdge, false, edge_tol);

  // Compute the fractional value of the closed curve
  const int n = approximating_polygon.numVertices();
  const double closure_wn =
    detail::linear_winding_number(q,
                                  approximating_polygon[n - 1],
                                  approximating_polygon[0],
                                  edge_tol);

  // If the point is on the boundary of the approximating polygon,
  //  or coincident with the curve (rare), then winding_number<polygon>
  //  doesn't return the right half-integer. Have to go edge-by-edge.
  if(isCoincident || isOnEdge)
  {
    closed_curve_wn = closure_wn;
    for(int i = 1; i < n; ++i)
    {
      closed_curve_wn +=
        detail::linear_winding_number(q,
                                      approximating_polygon[i - 1],
                                      approximating_polygon[i],
                                      edge_tol);
    }
  }

  return gwn + closed_curve_wn - closure_wn;
}

/*!
 * \brief Computes the GWN for a 2D point wrt a 2D NURBS curve
 *
 * \param [in] query The query point to test
 * \param [in] n The NURBS curve object 
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN by decomposing into rational Bezier curves
 *  and summing the resulting GWNs. Far-away curves can be evaluated
 *  without decomposition using direct formula.
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const NURBSCurve<T, 2>& n,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  const int deg = n.getDegree();
  if(deg <= 0) return 0.0;

  // Early return is possible for most points + curves
  if(!n.boundingBox().expand(edge_tol).contains(q))
  {
    return detail::linear_winding_number(q,
                                         n[0],
                                         n[n.getNumControlPoints() - 1],
                                         edge_tol);
  }

  // Decompose the NURBS curve into Bezier segments
  auto beziers = n.extractBezier();

  // Compute the GWN for each Bezier segment
  double gwn = 0.0;
  for(int i = 0; i < beziers.size(); i++)
  {
    gwn += winding_number(q, beziers[i], edge_tol, EPS);
  }

  return gwn;
}

/*!
 * \brief Computes the GWN for a 2D point wrt to a 2D curved polygon
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The CurvedPolygon object
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN for the curved polygon by summing the GWN for each curved edge
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const CurvedPolygon<T, 2>& cpoly,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  double ret_val = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    ret_val += winding_number(q, cpoly[i], edge_tol, EPS);
  }

  return ret_val;
}

/*!
 * \brief Computes the GWN for a 2D point wrt to a collection of 2D Bezier curves
 *
 * \param [in] query The query point to test
 * \param [in] carray The array of Bezier curves
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Sums the GWN at `query` for each curved edge
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const axom::Array<BezierCurve<T, 2>>& carray,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  double ret_val = 0.0;
  for(int i = 0; i < carray.size(); i++)
  {
    ret_val += winding_number(q, carray[i], edge_tol, EPS);
  }

  return ret_val;
}

/*!
 * \brief Computes the GWN for a 2D point wrt to a collection of 2D NURBS curves
 *
 * \param [in] query The query point to test
 * \param [in] narray The array of NURBS curves
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Sums the GWN at `query` for each curved edge
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const axom::Array<NURBSCurve<T, 2>>& narray,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  double ret_val = 0.0;
  for(int i = 0; i < narray.size(); i++)
  {
    ret_val += winding_number(q, narray[i], edge_tol, EPS);
  }

  return ret_val;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_H_
