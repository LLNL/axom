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
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/evaluate_integral.hpp"
#include "axom/primal/operators/detail/winding_number_impl.hpp"

// C++ includes
#include <cmath>
#include <stdio.h>
#include <iostream>

namespace axom
{
namespace primal
{

//@{
//! @name Winding number operations between 2D points and primatives

/*
 * \brief Compute the winding number with respect to a 2D line segment
 *
 * \param [in] q The query point to test
 * \param [in] s The line segment
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * \return double The generalized winding number
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const Segment<T, 2>& s,
                      double edge_tol = 1e-8)
{
  return detail::linear_winding_number(q, s[0], s[1], edge_tol);
}

/*
 * \brief Compute the winding number with respect to a 2D triangle
 *
 * \param [in] q The query point to test
 * \param [in] tri The triangle
 * \param [in] includeBoundary If true, points on the boundary are considered interior.
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * The triangle is assumed to be closed, so the winding number is an integer
 * 
 * \return int The integer winding number
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
 * \brief Computes the winding number for a point and a 2D polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [in] includeBoundary If true, points on the boundary are considered interior.
 * \param [in] EPS The tolerance level for collinearity
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
                   bool includeBoundary = false,
                   double EPS = 1e-8)
{
  const int nverts = P.numVertices();

  // If the query is a vertex, return a value interpreted
  //  as "inside" by evenodd or nonzero protocols
  if(axom::utilities::isNearlyEqual(P[0][0], R[0], EPS) &&
     axom::utilities::isNearlyEqual(P[0][1], R[1], EPS))
  {
    return includeBoundary;
  }

  int winding_num = 0;
  for(int i = 0; i < nverts; i++)
  {
    int j = (i == nverts - 1) ? 0 : i + 1;

    if(axom::utilities::isNearlyEqual(P[j][1], R[1], EPS))
    {
      if(axom::utilities::isNearlyEqual(P[j][0], R[0], EPS))
      {
        return includeBoundary;  // On vertex
      }
      else if(P[i][1] == R[1] && ((P[j][0] > R[0]) == (P[i][0] < R[0])))
      {
        return includeBoundary;  // On horizontal edge
      }
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

          // On edge
          if(axom::utilities::isNearlyEqual(det, 0.0, EPS))
          {
            return includeBoundary;
          }

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

          // On edge
          if(axom::utilities::isNearlyEqual(det, 0.0, EPS))
          {
            return includeBoundary;
          }

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
 * \brief Computes the generalized winding number for a single 2D Bezier curve
 *
 * \param [in] query The query point to test
 * \param [in] c The Bezier curve object 
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the winding number using a recursive, bisection algorithm,
 * using nearly-linear Bezier curves as a base case.
 * 
 * \return double the generalized winding number.
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
 * \brief Computes the generalized winding number for a 2D curved polygon
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The CurvedPolygon object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the winding number by summing the winding number for each curve
 * 
 * \return double the generalized winding number.
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
    ret_val +=
      detail::curve_winding_number_recursive(q, cpoly[i], false, edge_tol, EPS);
  }

  return ret_val;
}

//@}

//@{
//! @name Winding number operations between 3D points and primatives

/*!
 * \brief Computes the solid angle winding number for a 3D triangle
 *
 * \param [in] query The query point to test
 * \param [in] tri The 3D Triangle object
 * \param [in] isOnFace An optional return parameter if the point is on the triangle
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the winding number using the formula from 
 *  Oosterom, Strackee, "The Solid Angle of a Plane Triangle" 
 *  IEEE Transactions on Biomedical Engineering, Vol BME-30, No. 2, February 1983
 * with extra adjustments if the triangle takes up a full octant
 * 
 * \return double the generalized winding number.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Triangle<T, 3>& tri,
                      bool& isOnFace,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  if(tri.area() == 0) return 0;

  Vector<T, 3> a(q, tri[0]), b(q, tri[1]), c(q, tri[2]);

  // Compute norms. Possibly return early
  const double a_norm = a.norm();
  const double b_norm = b.norm();
  const double c_norm = c.norm();

  if(a_norm < edge_tol || b_norm < edge_tol || c_norm < edge_tol) return 0;

  double num = Vector<T, 3>::scalar_triple_product(a, b, c);
  if(axom::utilities::isNearlyEqual(num, 0.0, EPS))
  {
    isOnFace = true;
    return 0;
  }

  double denom = a_norm * b_norm * c_norm +
    a_norm * Vector<T, 3>::dot_product(b, c) +
    b_norm * Vector<T, 3>::dot_product(a, c) +
    c_norm * Vector<T, 3>::dot_product(a, b);

  // Handle direct cases where argument to atan is undefined
  if(axom::utilities::isNearlyEqual(denom, 0.0, EPS))
  {
    return (num > 0) ? 0.25 : -0.25;
  }

  if(denom > 0)
  {
    return 0.5 * M_1_PI * atan(num / denom);
  }
  else
  {
    if(num > 0)
    {
      return 0.5 * M_1_PI * atan(num / denom) + 0.5;
    }
    if(num < 0)
    {
      return 0.5 * M_1_PI * atan(num / denom) - 0.5;
    }
  }
}
/*!
 * \brief Computes the solid angle winding number for a 3D triangle
 *
 * Overload function without additional returning parameter
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Triangle<T, 3>& tri,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  bool isOnFace;
  return winding_number(q, tri, isOnFace, edge_tol, EPS);
}

/*!
 * \brief Computes the solid angle winding number for a 3D planar polygon
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
 * Triangulates the polygon and computes the triangular solid angle for each part
 * 
 * \return double the generalized winding number.
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
 * \brief Computes the solid angle winding number for a 3D planar polygon
 *
 * Overload function without additional returning parameter
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Polygon<T, 3>& poly,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  bool isOnFace;
  return winding_number(q, poly, isOnFace, edge_tol, EPS);
}

/*!
 * \brief Computes the solid angle winding number for a 3D convex polyhedron
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
 * Computes the faces of the polyhedron and computes the winding number for each.
 *
 * \return int The integer winding number.
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

template <typename T>
double winding_number_adapt_stokes(Point<T, 3> query,
                                   const BezierPatch<T>& bPatch,
                                   int& npts,
                                   bool includeBoundary = false,
                                   double quad_tol = 1e-8,
                                   double edge_tol = 1e-8,
                                   double EPS = 1e-8)
{
  const int quad_order = 20;
  const int ord_u = bPatch.getOrder_u();
  const int ord_v = bPatch.getOrder_v();
  const bool patchIsRational = bPatch.isRational();

  // Define three "axis-aligned" vector fields that
  //  integrate to the winding number
  auto x_field = [&query](Point<T, 3> x) -> Vector<T, 3> {
    Vector<T, 3> xmq = Vector<T, 3>(query, x);
    double denom1 = xmq[1] * xmq[1] + xmq[2] * xmq[2];
    double denom2 = xmq.norm();
    return Vector<T, 3>({0.0,
                         xmq[2] * xmq[0] / denom1 / denom2,
                         -xmq[1] * xmq[0] / denom1 / denom2});
  };

  auto y_field = [&query](Point<T, 3> x) -> Vector<T, 3> {
    Vector<T, 3> xmq = Vector<T, 3>(query, x);
    double denom1 = xmq[0] * xmq[0] + xmq[2] * xmq[2];
    double denom2 = xmq.norm();
    return Vector<T, 3>({-xmq[2] * xmq[1] / denom1 / denom2,
                         0.0,
                         xmq[0] * xmq[1] / denom1 / denom2});
  };

  auto z_field = [&query](Point<T, 3> x) -> Vector<T, 3> {
    Vector<T, 3> xmq = Vector<T, 3>(query, x);
    double denom1 = xmq[0] * xmq[0] + xmq[1] * xmq[1];
    double denom2 = xmq.norm();
    return Vector<T, 3>({xmq[1] * xmq[2] / denom1 / denom2,
                         -xmq[0] * xmq[2] / denom1 / denom2,
                         0.0});
  };

  // Define a bounding box that surrounds the patch
  BoundingBox<T, 3> bBox(bPatch.boundingBox());
  Point<T, 3> bBox_min = bBox.getMin();
  Point<T, 3> bBox_max = bBox.getMax();

  // Define a CurvedPolygon for the edges of the patch
  CurvedPolygon<T, 3> boundingPoly;
  boundingPoly.addEdge(bPatch.isocurve_u(0));
  boundingPoly.addEdge(bPatch.isocurve_v(1));
  boundingPoly.addEdge(bPatch.isocurve_u(1));
  boundingPoly.addEdge(bPatch.isocurve_v(0));
  boundingPoly[0].reverseOrientation();
  boundingPoly[1].reverseOrientation();

  // Check if the (z-q)-axis is inside the bounding box
  if(!(bBox_min[0] < query[0] && query[0] < bBox_max[0] &&
       bBox_min[1] < query[1] && query[1] < bBox_max[1]))
  {
    return 0.25 * M_1_PI *
      evaluate_vector_line_integral(boundingPoly, z_field, quad_order, npts);
  }

  // Check if the (y-q)-axis is inside the bounding box
  if(!(bBox_min[0] < query[0] && query[0] < bBox_max[0] &&
       bBox_min[2] < query[2] && query[2] < bBox_max[2]))
  {
    return 0.25 * M_1_PI *
      evaluate_vector_line_integral(boundingPoly, y_field, quad_order, npts);
  }

  // Check if the (x-q)-axis is inside the bounding box
  if(!(bBox_min[1] < query[1] && query[1] < bBox_max[1] &&
       bBox_min[2] < query[2] && query[2] < bBox_max[2]))
  {
    return 0.25 * M_1_PI *
      evaluate_vector_line_integral(boundingPoly, x_field, quad_order, npts);
  }

  // If all of these tests fall through, then we are inside the bounding box
  //  and need to split
  BezierPatch<T> p1, p2, p3, p4;
  bPatch.split(0.5, 0.5, p1, p2, p3, p4);
  // clang-format off
  return winding_number_adapt_stokes(query, p1, npts, includeBoundary, quad_tol, edge_tol, EPS) +
    winding_number_adapt_stokes(query, p2, npts, includeBoundary, quad_tol, edge_tol, EPS) +
    winding_number_adapt_stokes(query, p3, npts, includeBoundary, quad_tol, edge_tol, EPS) +
    winding_number_adapt_stokes(query, p4, npts, includeBoundary, quad_tol, edge_tol, EPS);
  // clang-format on

  //if(axom::utilities::isNearlyEqual(wn + closure_wn, 0.0, quad_tol))
  //{
  //  // bPatch.python_print(std::cout, npts, true);
  //  //std::cout << wn << std::endl;
  //  return wn;
  //}
  //else  // If not, split recursively until it is
  //{
  //  BezierPatch<T> p1, p2, p3, p4;
  //  bPatch.split(0.5, 0.5, p1, p2, p3, p4);
  //  return winding_number(query, p1, npts, includeBoundary, quad_tol, edge_tol, EPS) +
  //    winding_number(query, p2, npts, includeBoundary, quad_tol, edge_tol, EPS) +
  //    winding_number(query, p3, npts, includeBoundary, quad_tol, edge_tol, EPS) +
  //    winding_number(query, p4, npts, includeBoundary, quad_tol, edge_tol, EPS);
  //}
}

template <typename T>
double winding_number_mfem_triangle(Point<T, 3> query,
                                    const BezierPatch<T>& bPatch,
                                    int& npts,
                                    bool includeBoundary = false,
                                    double quad_tol = 1e-8,
                                    double edge_tol = 1e-8,
                                    double EPS = 1e-8)
{
  const int quad_order = 20;
  const int ord_u = bPatch.getOrder_u();
  const int ord_v = bPatch.getOrder_v();
  const bool patchIsRational = bPatch.isRational();

  BoundingBox<T, 3> bBox(bPatch.boundingBox());
  if(bBox.contains(query))
  {
    BezierPatch<T> p1, p2, p3, p4;
    bPatch.split(0.5, 0.5, p1, p2, p3, p4);
    // clang-format off
    return winding_number_mfem_triangle(query, p1, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_mfem_triangle(query, p2, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_mfem_triangle(query, p3, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_mfem_triangle(query, p4, npts, includeBoundary, quad_tol, edge_tol, EPS);
    // clang-format on
  }

  // Define an integrand for computing winding numbers of square patches
  BezierPatch<T> surface_patch;
  auto wn_integrand = [&query, &surface_patch](double u, double v) -> double {
    Vector<T, 3> xmq(query, surface_patch.evaluate(u, v));
    return 0.25 * M_1_PI *
      Vector<T, 3>::dot_product(xmq, surface_patch.normal(u, v)) /
      std::pow(xmq.norm(), 3);
  };

  // Apply the integrand to the main patch
  surface_patch = bPatch;
  double wn =
    evaluate_parameter_integral_square(wn_integrand, 2 * quad_order, npts);

  // Define a pyramid-like object whose top faces are a ruled surface
  //  to the edges of the patch
  Point<T, 3> cone_vertex = bBox.getCentroid();

  BezierCurve<T, 3> edge_curves[4];
  BezierCurve<T, 3> edge_curve;
  auto wn_integrand_tri =
    [&query, &edge_curve, &cone_vertex](double u, double v) -> double {
    double t = u / (u + v);
    Point<T, 3> c = edge_curve.evaluate(t);
    Point<T, 3> node = cone_vertex + (u + v) * (c - cone_vertex);

    Vector<T, 3> c_prime = edge_curve.dt(t);
    Vector<T, 3> tangent_u = v * c_prime / (u + v) + c - cone_vertex;
    Vector<T, 3> tangent_v = -u * c_prime / (u + v) + c - cone_vertex;
    Vector<T, 3> normal = Vector<T, 3>::cross_product(tangent_u, tangent_v);

    Vector<T, 3> xmq(query, node);
    return 0.25 * M_1_PI * Vector<T, 3>::dot_product(xmq, normal) /
      std::pow(xmq.norm(), 3);
  };

  // Create ruled surfaces for the (0, v) and (1, v) isocurves
  edge_curves[0].setOrder(ord_v);
  edge_curves[2].setOrder(ord_v);
  if(patchIsRational)
  {
    edge_curves[0].makeRational();
    edge_curves[2].makeRational();
  }
  for(int i = 0; i <= ord_v; ++i)
  {
    edge_curves[0][i] = bPatch(ord_v, i);
    edge_curves[2][i] = bPatch(0, ord_v - i);
    if(patchIsRational)
    {
      edge_curves[0].setWeight(i, bPatch.getWeight(ord_v, i));
      edge_curves[2].setWeight(i, bPatch.getWeight(0, ord_v - i));
    }
  }

  // Create ruled surfaces for the (u, 0) and (u, 1) isocurves
  edge_curves[1].setOrder(ord_u);
  edge_curves[3].setOrder(ord_u);
  if(patchIsRational)
  {
    edge_curves[1].makeRational();
    edge_curves[3].makeRational();
  }
  for(int i = 0; i <= ord_u; ++i)
  {
    edge_curves[1][i] = bPatch(ord_u - i, ord_u);
    edge_curves[3][i] = bPatch(i, 0);
    if(patchIsRational)
    {
      edge_curves[1].setWeight(i, bPatch.getWeight(ord_u - i, ord_u));
      edge_curves[3].setWeight(i, bPatch.getWeight(i, 0));
    }
  }

  // Apply the special triangular quadrature to each face
  double closure_wn = 0;
  for(int i = 0; i < 4; ++i)
  {
    edge_curve = edge_curves[i];
    closure_wn +=
      evaluate_parameter_integral_triangle(wn_integrand_tri, quad_order, npts);
  }

  if(axom::utilities::isNearlyEqual(wn + closure_wn, 0.0, quad_tol))
  {
    return wn;
  }
  else  // If not, split recursively until it is
  {
    BezierPatch<T> p1, p2, p3, p4;
    bPatch.split(0.5, 0.5, p1, p2, p3, p4);
    // clang-format off
    return winding_number_mfem_triangle(query, p1, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_mfem_triangle(query, p2, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_mfem_triangle(query, p3, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_mfem_triangle(query, p4, npts, includeBoundary, quad_tol, edge_tol, EPS);
    // clang-format on
  }
}

template <typename T>
double winding_number_gauss_tensor(Point<T, 3> query,
                                   const BezierPatch<T>& bPatch,
                                   int& npts,
                                   bool includeBoundary = false,
                                   double quad_tol = 1e-8,
                                   double edge_tol = 1e-8,
                                   double EPS = 1e-8)
{
  const int quad_order = 9;
  const int ord_u = bPatch.getOrder_u();
  const int ord_v = bPatch.getOrder_v();
  const bool patchIsRational = bPatch.isRational();

  BoundingBox<T, 3> bBox(bPatch.boundingBox());
  if(bBox.contains(query))
  {
    BezierPatch<T> p1, p2, p3, p4;
    bPatch.split(0.5, 0.5, p1, p2, p3, p4);
    // clang-format off
    return winding_number_gauss_tensor(query, p1, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_gauss_tensor(query, p2, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_gauss_tensor(query, p3, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_gauss_tensor(query, p4, npts, includeBoundary, quad_tol, edge_tol, EPS);
    // clang-format on
  }

  // Define an integrand for computing winding numbers of square patches
  BezierPatch<T> surface_patch;
  auto wn_integrand = [&query, &surface_patch](double u, double v) -> double {
    Vector<T, 3> xmq(query, surface_patch.evaluate(u, v));
    return 0.25 * M_1_PI *
      Vector<T, 3>::dot_product(xmq, surface_patch.normal(u, v)) /
      std::pow(xmq.norm(), 3);
  };

  // Apply the integrand to the main patch
  surface_patch = bPatch;
  double wn =
    evaluate_parameter_integral_square(wn_integrand, 2 * quad_order, npts);

  return wn;

  // Define a pyramid-like object whose top faces are a ruled surface
  //  to the edges of the patch
  Point<T, 3> cone_vertex = bBox.getCentroid();

  /* ==================== Use Patch algorithm for the cone ==================== */
  BezierPatch<T> cone_faces[4];

  // Create degenerate patch for the (0, v) and (1, v) isocurves
  cone_faces[0].setOrder(ord_v, 1);
  cone_faces[2].setOrder(ord_v, 1);
  if(patchIsRational)
  {
    cone_faces[0].makeRational();
    cone_faces[2].makeRational();
  }
  for(int i = 0; i <= ord_v; ++i)
  {
    cone_faces[0](i, 0) = cone_vertex;
    cone_faces[0](i, 1) = bPatch(ord_v, i);

    cone_faces[2](i, 0) = cone_vertex;
    cone_faces[2](i, 1) = bPatch(0, ord_v - i);

    if(patchIsRational)
    {
      cone_faces[0].setWeight(i, 1, bPatch.getWeight(ord_v, i));
      cone_faces[2].setWeight(i, 1, bPatch.getWeight(0, ord_v - i));
    }
  }

  // Create degenerate patch for the (u, 0) and (u, 1) isocurves
  cone_faces[1].setOrder(ord_u, 1);
  cone_faces[3].setOrder(ord_u, 1);
  if(patchIsRational)
  {
    cone_faces[1].makeRational();
    cone_faces[3].makeRational();
  }
  for(int i = 0; i <= ord_u; ++i)
  {
    cone_faces[1](i, 0) = cone_vertex;
    cone_faces[1](i, 1) = bPatch(ord_u - i, ord_u);

    cone_faces[3](i, 0) = cone_vertex;
    cone_faces[3](i, 1) = bPatch(i, 0);

    if(patchIsRational)
    {
      cone_faces[1].setWeight(i, 1, bPatch.getWeight(ord_u - i, ord_u));
      cone_faces[3].setWeight(i, 1, bPatch.getWeight(i, 0));
    }
  }

  // Apply the integrand to the four cone faces
  double closure_wn = 0;
  for(int n = 0; n < 4; ++n)
  {
    surface_patch = cone_faces[n];
    closure_wn +=
      evaluate_parameter_integral_square(wn_integrand, quad_order, npts);
  }

  if(axom::utilities::isNearlyEqual(wn + closure_wn, 0.0, quad_tol))
  {
    return wn;
  }
  else  // If not, split recursively until it is
  {
    BezierPatch<T> p1, p2, p3, p4;
    bPatch.split(0.5, 0.5, p1, p2, p3, p4);
    // clang-format off
    return winding_number_gauss_tensor(query, p1, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_gauss_tensor(query, p2, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_gauss_tensor(query, p3, npts, includeBoundary, quad_tol, edge_tol, EPS) +
      winding_number_gauss_tensor(query, p4, npts, includeBoundary, quad_tol, edge_tol, EPS);
    // clang-format on
  }
}

//@}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_H_
