// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file winding_number.hpp
 *
 * \brief Consists of methods to compute the generalized winding number (GWN) 
 *        for points with respect to various geometric objects.
 */

#ifndef AXOM_PRIMAL_WINDING_NUMBER_HPP_
#define AXOM_PRIMAL_WINDING_NUMBER_HPP_

// Axom includes
#include "axom/core.hpp"
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/detail/winding_number_impl.hpp"

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
//@{
//! @name Winding number operations between 2D points and primitives

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
 * return the GWN of the original curve. (See )
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

  // Early return is possible for must points + curves
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
 * \param [in] cpoly The CurvedPolygon object
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
  AXOM_UNUSED_VAR(EPS);

  double ret_val = 0.0;
  for(int i = 0; i < carray.size(); i++)
  {
    ret_val += winding_number(q, carray[i], false, edge_tol);
  }

  return ret_val;
}
//@}

//@{
//! @name Winding number operations between 3D points and primitives

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D triangle
 *
 * \param [in] query The query point to test
 * \param [in] tri The 3D Triangle object
 * \param [in] isOnFace An optional return parameter if the point is on the triangle
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN as the solid angle modulo 4pi using the formula from 
 *  Oosterom, Strackee, "The Solid Angle of a Plane Triangle" 
 *  IEEE Transactions on Biomedical Engineering, Vol BME-30, No. 2, February 1983
 * with extra adjustments if the triangle takes up a full octant
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Triangle<T, 3>& tri,
                      bool& isOnFace,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  using Vec3 = Vector<T, 3>;

  if(tri.area() == 0)
  {
    return 0;
  }

  const Vec3 a = tri[0] - q;
  const Vec3 b = tri[1] - q;
  const Vec3 c = tri[2] - q;

  // Compute norms. Possibly return early
  const double a_norm = a.norm();
  const double b_norm = b.norm();
  const double c_norm = c.norm();

  if(a_norm < edge_tol || b_norm < edge_tol || c_norm < edge_tol)
  {
    return 0;
  }

  const double num = Vec3::scalar_triple_product(a, b, c);
  if(axom::utilities::isNearlyEqual(num, 0.0, EPS))
  {
    isOnFace = true;
    return 0;
  }

  const double denom = a_norm * b_norm * c_norm + a_norm * b.dot(c) +
    b_norm * a.dot(c) + c_norm * a.dot(b);

  // Handle direct cases where argument to atan is undefined
  if(axom::utilities::isNearlyEqual(denom, 0.0, EPS))
  {
    return (num > 0) ? 0.25 : -0.25;
  }

  // Note: denom==0 and num==0 handled above
  if(denom > 0)
  {
    return 0.5 * M_1_PI * atan(num / denom);
  }
  else
  {
    return (num > 0) ? 0.5 * M_1_PI * atan(num / denom) + 0.5
                     : 0.5 * M_1_PI * atan(num / denom) - 0.5;
  }
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D triangle
 *
 * \param [in] query The query point to test
 * \param [in] tri The 3D Triangle object
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN for the triangle without an additional return parameter
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Triangle<T, 3>& tri,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  bool isOnFace = false;
  return winding_number(q, tri, isOnFace, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D planar polygon
 *
 * \param [in] query The query point to test
 * \param [in] poly The Polygon object
 * \param [in] isOnFace Return variable to show if the point is on the polygon
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \pre Assumes the polygon is planar. Otherwise, a meaningless value is returned.
 * 
 * Triangulates the polygon and computes the triangular GWN for each component
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Polygon<T, 3>& poly,
                      bool& isOnFace,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  const int num_verts = poly.numVertices();
  if(num_verts < 3)
  {
    return 0;
  }

  double wn = 0.0;
  for(int i = 0; i < num_verts - 2; ++i)
  {
    wn += winding_number(q,
                         Triangle<T, 3>(poly[0], poly[i + 1], poly[i + 2]),
                         isOnFace,
                         edge_tol,
                         EPS);
  }

  return wn;
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D planar polygon
 *
 * \param [in] query The query point to test
 * \param [in] poly The Polygon object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \pre Assumes the polygon is planar. Otherwise, a meaningless value is returned.
 * 
 * Computes the GWN for the polygon without an additional return parameter
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Polygon<T, 3>& poly,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  bool isOnFace = false;
  return winding_number(q, poly, isOnFace, edge_tol, EPS);
}

/*!
 * \brief Computes the winding number for a 3D point wrt a 3D convex polyhedron
 *
 * \param [in] query The query point to test
 * \param [in] poly The Polyhedron object
 * \param [in] includeBoundary If true, points on the boundary are considered interior.
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \pre Expects the polyhedron to be convex and closed so that the returned value is an integer.
 * 
 * Computes the faces of the polyhedron and computes the GWN for each.
 * The sum is then rounded to the nearest integer, as the shape is assumed to be closed.
 * 
 * \return The integer winding number.
 */
template <typename T>
int winding_number(const Point<T, 3>& query,
                   const Polyhedron<T, 3>& poly,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8,
                   double EPS = 1e-8)
{
  SLIC_ASSERT(poly.hasNeighbors());
  const int num_verts = poly.numVertices();

  axom::Array<int> faces(num_verts * num_verts), face_size(2 * num_verts),
    face_offset(2 * num_verts);
  int face_count;

  poly.getFaces(faces.data(), face_size.data(), face_offset.data(), face_count);

  bool isOnFace = false;
  double wn = 0;
  for(int i = 0; i < face_count; ++i)
  {
    const int N = face_size[i];
    const int i_offset = face_offset[i];
    Polygon<T, 3> the_face(N);
    for(int j = 0; j < N; ++j)
    {
      the_face.addVertex(poly[faces[i_offset + j]]);
    }

    wn += winding_number(query, the_face, isOnFace, edge_tol, EPS);

    if(isOnFace)
    {
      return includeBoundary;
    }
  }

  return std::lround(wn);
}

#ifdef AXOM_USE_MFEM

/*
 * \brief Computes the GWN for a 3D point wrt a 3D Bezier patch
 *
 * \param [in] query The query point to test
 * \param [in] bPatch The Bezier patch object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] quad_tol The maximum relative error allowed in the quadrature
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * \param [in] depth The current recursive depth
 * 
 * Computes the generalized winding number for a Bezier patch using Stokes theorem.
 *
 * \note Warning: This algorithm is only tested to high accuracy for queries within
 *  1e-5 of the surface. Otherwise, it will return less accurate results.
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& query,
                      const BezierPatch<T, 3>& bPatch,
                      const double edge_tol = 1e-8,
                      const double quad_tol = 1e-8,
                      const double EPS = 1e-8,
                      const int depth = 0)
{
  const int ord_u = bPatch.getOrder_u();
  const int ord_v = bPatch.getOrder_v();
  const bool patchIsRational = bPatch.isRational();
  const double edge_tol_sq = edge_tol * edge_tol;

  // Fix the number of quadrature nodes arbitrarily, but high enough
  //  to `catch` near singularities for refinement
  constexpr int quad_npts = 30;

  // Early return if the patch is approximately polygonal.
  //  Very slight variations in curvature requires small EPS tolerance
  constexpr int MAX_DEPTH = 10;
  if(depth >= MAX_DEPTH || bPatch.isPolygonal(EPS))
  {
    return winding_number(
      query,
      Polygon<T, 3>(axom::Array<Point<T, 3>>(
        {bPatch(0, 0), bPatch(ord_u, 0), bPatch(ord_u, ord_v), bPatch(0, ord_v)})),
      edge_tol,
      PRIMAL_TINY);
  }

  // Use a specific kind of recursion if we are within tol of an endpoint.
  //  Split the surface closer to the corner, assume smallest patch is polygonal,
  //  and set a new edge_tol so corners of the new patch aren't marked as coincident
  constexpr double edge_offset = 0.01;
  if(squared_distance(query, bPatch(0, 0)) <= edge_tol_sq)
  {
    BezierPatch<T, 3> p1, p2, p3, p4;
    bPatch.split(0.0 + edge_offset, 0.0 + edge_offset, p1, p2, p3, p4);
    double new_edge_tol = 0.5 *
      sqrt(axom::utilities::min(
        squared_distance(query, bPatch.evaluate(0.0, 0.0 + edge_offset)),
        squared_distance(query, bPatch.evaluate(0.0 + edge_offset, 0.0))));
    new_edge_tol = axom::utilities::min(new_edge_tol, edge_tol);

    return winding_number(query, p2, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p3, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p4, new_edge_tol, quad_tol, EPS, depth + 1);
  }
  if(squared_distance(query, bPatch(ord_u, 0)) <= edge_tol_sq)
  {
    BezierPatch<T, 3> p1, p2, p3, p4;
    bPatch.split(1.0 - edge_offset, 0.0 + edge_offset, p1, p2, p3, p4);
    double new_edge_tol = 0.5 *
      sqrt(axom::utilities::min(
        squared_distance(query, bPatch.evaluate(1.0, 0.0 + edge_offset)),
        squared_distance(query, bPatch.evaluate(1.0 - edge_offset, 0.0))));
    new_edge_tol = axom::utilities::min(new_edge_tol, edge_tol);

    return winding_number(query, p1, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p3, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p4, new_edge_tol, quad_tol, EPS, depth + 1);
  }
  if(squared_distance(query, bPatch(0, ord_v)) <= edge_tol_sq)
  {
    BezierPatch<T, 3> p1, p2, p3, p4;
    bPatch.split(0.0 + edge_offset, 1.0 - edge_offset, p1, p2, p3, p4);
    double new_edge_tol = 0.5 *
      sqrt(axom::utilities::min(
        squared_distance(query, bPatch.evaluate(0.0 + edge_offset, 1.0)),
        squared_distance(query, bPatch.evaluate(0.0, 1.0 - edge_offset))));
    new_edge_tol = axom::utilities::min(new_edge_tol, edge_tol);

    return winding_number(query, p1, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p2, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p4, new_edge_tol, quad_tol, EPS, depth + 1);
  }
  if(squared_distance(query, bPatch(ord_u, ord_v)) <= edge_tol_sq)
  {
    BezierPatch<T, 3> p1, p2, p3, p4;
    bPatch.split(1.0 - edge_offset, 1.0 - edge_offset, p1, p2, p3, p4);
    double new_edge_tol = 0.5 *
      sqrt(axom::utilities::min(
        squared_distance(query, bPatch.evaluate(1.0, 1.0 - edge_offset)),
        squared_distance(query, bPatch.evaluate(1.0 - edge_offset, 1.0))));
    new_edge_tol = axom::utilities::min(new_edge_tol, edge_tol);

    return winding_number(query, p1, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p2, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p3, new_edge_tol, quad_tol, EPS, depth + 1);
  }

  /* 
   * To use Stokes theorem, we need to identify a separating plane between
   * `query` and the surface, guaranteed through a bounding box.
   * If it does, need to do geometric refinement: Splitting and rotating the curve
   * until we can guarantee this.
   */
  CurvedPolygon<T, 3> boundingPoly(4);

  // Define vector fields whose curl gives us the winding number
  detail::SingularityAxis field_direction;

  // Check an axis-aligned bounding box (most surfaces satisfy this condition)
  BoundingBox<T, 3> bBox(bPatch.boundingBox().expand(edge_tol));
  const bool exterior_x =
    bBox.getMin()[0] > query[0] || query[0] > bBox.getMax()[0];
  const bool exterior_y =
    bBox.getMin()[1] > query[1] || query[1] > bBox.getMax()[1];
  const bool exterior_z =
    bBox.getMin()[2] > query[2] || query[2] > bBox.getMax()[2];

  if(exterior_y || exterior_z)
  {
    field_direction = detail::SingularityAxis::x;
  }
  else if(exterior_x || exterior_z)
  {
    field_direction = detail::SingularityAxis::y;
  }
  else if(exterior_x || exterior_y)
  {
    field_direction = detail::SingularityAxis::z;
  }
  else
  {
    // Next, check an oriented bounding box.
    // If we are interior to the oriented bounding box, then we
    //  cannot guarantee a separating plane, and need geometric refinement.
    OrientedBoundingBox<T, 3> oBox(bPatch.orientedBoundingBox().expand(edge_tol));
    if(oBox.contains(query))
    {
      BezierPatch<T, 3> p1, p2, p3, p4;
      bPatch.split(0.5, 0.5, p1, p2, p3, p4);
      return winding_number(query, p1, edge_tol, quad_tol, EPS, depth + 1) +
        winding_number(query, p2, edge_tol, quad_tol, EPS, depth + 1) +
        winding_number(query, p3, edge_tol, quad_tol, EPS, depth + 1) +
        winding_number(query, p4, edge_tol, quad_tol, EPS, depth + 1);
    }

    // Otherwise, we can apply a rotation to a z-aligned field.
    field_direction = detail::SingularityAxis::rotated;

    // Lambda to generate a 3D rotation matrix from an angle and axis
    // Formulation from https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    auto angleAxisRotMatrix =
      [](double theta, const Vector<T, 3>& axis) -> numerics::Matrix<T> {
      const auto unitized = axis.unitVector();
      const double x = unitized[0], y = unitized[1], z = unitized[2];
      const double c = cos(theta), s = sin(theta), C = 1 - c;

      auto matx = numerics::Matrix<T>::zeros(3, 3);

      matx(0, 0) = x * x * C + c;
      matx(0, 1) = x * y * C - z * s;
      matx(0, 2) = x * z * C + y * s;

      matx(1, 0) = y * x * C + z * s;
      matx(1, 1) = y * y * C + c;
      matx(1, 2) = y * z * C - x * s;

      matx(2, 0) = z * x * C - y * s;
      matx(2, 1) = z * y * C + x * s;
      matx(2, 2) = z * z * C + c;

      return matx;
    };

    // Lambda to rotate the input point using the provided rotation matrix
    auto rotate_point = [&query](const numerics::Matrix<T>& matx,
                                 const Point<T, 3> input) -> Point<T, 3> {
      Vector<T, 3> shifted(query, input);
      Vector<T, 3> rotated;
      numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
      return Point<T, 3>(
        {rotated[0] + query[0], rotated[1] + query[1], rotated[2] + query[2]});
    };

    // Find vector from query to the bounding box
    Point<T, 3> closest = closest_point(query, oBox);
    Vector<T, 3> v0 = Vector<T, 3>(query, closest).unitVector();

    // Find the direction of a ray perpendicular to that
    Vector<T, 3> v1;
    if(axom::utilities::isNearlyEqual(v0[0], v0[1], EPS))
    {
      v1 = Vector<T, 3>({v0[2], v0[2], -v0[0] - v0[1]}).unitVector();
    }
    else
    {
      v1 = Vector<T, 3>({-v0[1] - v0[2], v0[0], v0[0]}).unitVector();
    }

    // Rotate v0 around v1 until it is perpendicular to the plane spanned by k and v1
    double ang = (v0[2] < 0 ? 1.0 : -1.0) *
      acos(axom::utilities::clampVal(
        -(v0[0] * v1[1] - v0[1] * v1[0]) / sqrt(v1[0] * v1[0] + v1[1] * v1[1]),
        -1.0,
        1.0));
    auto rotator = angleAxisRotMatrix(ang, v1);

    // Collect rotated curves into the curved Polygon
    // Set up the (0, v) and (1, v) isocurves, rotated
    boundingPoly[0].setOrder(ord_v);
    boundingPoly[2].setOrder(ord_v);
    if(patchIsRational)
    {
      boundingPoly[0].makeRational();
      boundingPoly[2].makeRational();
    }
    for(int q = 0; q <= ord_v; ++q)
    {
      boundingPoly[0][q] = rotate_point(rotator, bPatch(ord_v, q));
      boundingPoly[2][q] = rotate_point(rotator, bPatch(0, ord_v - q));

      if(patchIsRational)
      {
        boundingPoly[0].setWeight(q, bPatch.getWeight(ord_v, q));
        boundingPoly[2].setWeight(q, bPatch.getWeight(0, ord_v - q));
      }
    }

    // Set up the (u, 0) and (u, 1) isocurves
    boundingPoly[1].setOrder(ord_u);
    boundingPoly[3].setOrder(ord_u);
    if(patchIsRational)
    {
      boundingPoly[1].makeRational();
      boundingPoly[3].makeRational();
    }
    for(int p = 0; p <= ord_u; ++p)
    {
      boundingPoly[1][p] = rotate_point(rotator, bPatch(ord_u - p, ord_u));
      boundingPoly[3][p] = rotate_point(rotator, bPatch(p, 0));

      if(patchIsRational)
      {
        boundingPoly[1].setWeight(p, bPatch.getWeight(ord_u - p, ord_u));
        boundingPoly[3].setWeight(p, bPatch.getWeight(p, 0));
      }
    }
  }

  // Set up the polygon if we don't need to do any rotation or splitting.
  if(field_direction != detail::SingularityAxis::rotated)
  {
    //  Add the relevant bounding curves to the patch.
    boundingPoly[0] = bPatch.isocurve_u(0);
    boundingPoly[0].reverseOrientation();

    boundingPoly[1] = bPatch.isocurve_v(1);
    boundingPoly[1].reverseOrientation();

    boundingPoly[2] = bPatch.isocurve_u(1);
    boundingPoly[3] = bPatch.isocurve_v(0);
  }

  // Iterate over the edges of the bounding curved polygon, add up the results
  double wn = 0;
  for(int n = 0; n < 4; ++n)
  {
    wn += detail::stokes_winding_number(query,
                                        boundingPoly[n],
                                        field_direction,
                                        quad_npts,
                                        quad_tol);
  }

  return wn;
}
#endif

//@}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_H_
