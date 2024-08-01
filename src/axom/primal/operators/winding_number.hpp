// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
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
                   bool& isOnEdge,
                   bool includeBoundary = false,
                   double EPS = 1e-8)
{
  const int nverts = P.numVertices();
  isOnEdge = false;

  // If the query is a vertex, return a value interpreted
  //  as "inside" by evenodd or nonzero protocols
  if(axom::utilities::isNearlyEqual(P[0][0], R[0], EPS) &&
     axom::utilities::isNearlyEqual(P[0][1], R[1], EPS))
  {
    isOnEdge = true;
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
        isOnEdge = true;
        return includeBoundary;  // On vertex
      }
      else if(P[i][1] == R[1] && ((P[j][0] > R[0]) == (P[i][0] < R[0])))
      {
        isOnEdge = true;
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
            isOnEdge = true;
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
            isOnEdge = true;
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
 * \brief Computes the solid angle winding number for a 2D polygon
 *
 * Overload function without additional returning parameter
 */
template <typename T>
int winding_number(const Point<T, 2>& R,
                   const Polygon<T, 2>& P,
                   bool includeBoundary = false,
                   double EPS = 1e-8)
{
  bool isOnEdge = false;
  return winding_number(R, P, isOnEdge, includeBoundary, EPS);
}

/*!
 * \brief Computes the generalized winding number for a single 2D Bezier curve
 *
 * \param [in] query The query point to test
 * \param [in] c The Bezier curve object 
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
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
                      int& nevals,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  return detail::curve_winding_number_recursive(q, c, false, nevals, edge_tol);
  const int ord = c.getOrder();
  if(ord <= 0) return 0.0;
  Polygon<T, 2> approxogon(2);
  approxogon.addVertex(c[0]);

  int dummy_int = 0;

  double wn = 0.0;
  //std::cout << "------------" << std::endl;
  detail::winding_number_adaptive_linear(q,
                                         c,
                                         false,
                                         dummy_int,
                                         edge_tol,
                                         edge_tol,
                                         approxogon,
                                         wn);
  //std::cout << "------------" << std::endl;
  approxogon.addVertex(c[ord]);

  //std::cout << approxogon.numVertices() << std::endl;

  return wn + detail::approxogon_winding_number(q, approxogon, edge_tol);
}

template <typename T>
double winding_number_quad(const Point<T, 2>& q,
                           const BezierCurve<T, 2>& c,
                           int npts)
{
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

  double quad_result = 0.0;
  for(int k = 0; k < quad.GetNPoints(); ++k)
  {
    Point<T, 2> x_quad = c.evaluate(quad.IntPoint(k).x);
    Vector<T, 2> dx_quad = c.dt(quad.IntPoint(k).x);

    double query_dist = squared_distance(x_quad, q);
    // clang-format off
    double query_orient = axom::numerics::determinant(x_quad[0] - q[0], dx_quad[0],
                                                      x_quad[1] - q[1], dx_quad[1]);
    // clang-format on

    quad_result += quad.IntPoint(k).weight * query_orient / query_dist;
  }

  return 0.5 * M_1_PI * quad_result;
}

template <typename T>
double winding_number_clipping(const Point<T, 2>& q,
                               const BezierCurve<T, 2>& c,
                               int& nevals,
                               double edge_tol = 1e-8,
                               double EPS = 1e-8)
{
  const int ord = c.getOrder();
  Segment<double, 2> closure(c[ord], c[0]);
  int inter = 0;

  BoundingBox<T, 2> bb = c.boundingBox();
  if(bb.contains(q))
  {
    // Check the 4 edges of the bounding box to find the closest
    double edge_dist = bb.getMax()[0] - q[0];
    Vector<double, 2> direction {1.0, 0.0};

    if(bb.getMax()[1] - q[1] < edge_dist)
    {
      direction = Vector<T, 2> {0.0, 1.0};
      edge_dist = bb.getMax()[1] - q[1];
    }
    if(q[0] - bb.getMin()[0] < edge_dist)
    {
      direction = Vector<T, 2> {-1.0, 0.0};
      edge_dist = q[0] - bb.getMin()[0];
    }
    if(q[1] - bb.getMin()[1] < edge_dist)
    {
      direction = Vector<T, 2> {0.0, -1.0};
    }

    Ray<T, 2> ray(q, direction);
    inter = detail::ray_casting_bezier_clipping(c, ray, nevals, edge_tol);
    double u, tmp;
    if(intersect(ray, closure, u, tmp, 0.0))
    {
      Vector<T, 2> e1 = ray.direction().unitVector();
      Vector<T, 2> e2 = Vector<T, 2> {-e1[1], e1[0]};

      double alpha = Vector<T, 2>(ray.origin(), c[ord]).dot(e2);

      (alpha < 0) ? inter++ : inter--;
    }
  }
  // If not contained in the bounding box, just do a return

  return inter - winding_number(q, closure, edge_tol);
}

template <typename T>
double winding_number_bisection(const Point<T, 2>& q,
                                const BezierCurve<T, 2>& c,
                                int& nevals,
                                double edge_tol = 1e-8,
                                double EPS = 1e-8)
{
  const int ord = c.getOrder();
  Segment<double, 2> closure(c[ord], c[0]);
  int crossing_num = 0;
  int inter = 0;

  BoundingBox<T, 2> bb = c.boundingBox();
  if(bb.contains(q))
  {
    // Check the 4 edges of the bounding box to find the closest
    double edge_dist = bb.getMax()[0] - q[0];
    Vector<double, 2> direction {1.0, 0.0};

    if(bb.getMax()[1] - q[1] < edge_dist)
    {
      direction = Vector<T, 2> {0.0, 1.0};
      edge_dist = bb.getMax()[1] - q[1];
    }
    if(q[0] - bb.getMin()[0] < edge_dist)
    {
      direction = Vector<T, 2> {-1.0, 0.0};
      edge_dist = q[0] - bb.getMin()[0];
    }
    if(q[1] - bb.getMin()[1] < edge_dist)
    {
      direction = Vector<T, 2> {0.0, -1.0};
    }

    Ray<T, 2> ray(q, direction);
    std::vector<double> cp;
    std::vector<double> rp;
    intersect(c, ray, cp, rp, nevals, edge_tol);

    for(auto c0 : cp)
    {
      Vector<T, 2> e1 = ray.direction().unitVector();
      Vector<T, 2> e2 = c.dt(c0);

      double determinant =
        axom::numerics::determinant(e1[0], e2[0], e1[1], e2[1]);

      (determinant > 0) ? crossing_num++ : crossing_num--;
    }

    double u, tmp;
    if(intersect(ray, closure, u, tmp, edge_tol))
    {
      Vector<T, 2> e1 = ray.direction().unitVector();
      Vector<T, 2> e2 = Vector<T, 2> {-e1[1], e1[0]};

      double alpha = Vector<T, 2>(q, c[ord]).dot(e2);

      (alpha < 0) ? crossing_num++ : crossing_num--;
    }
  }
  // If not contained in the bounding box, just do a return

  return crossing_num - winding_number(q, closure, edge_tol);
}

/// Compute the number of intersections that a ray makes with an array of Bezier curves
template <typename T>
int crossing_number(const Point<T, 2>& q,
                    const axom::Array<BezierCurve<T, 2>>& cs,
                    const BoundingBox<T, 2>& bb,
                    double edge_tol = 1e-8,
                    double EPS = 1e-8)
{
  int nevals = 0;
  int inter = 0;

  double edge_dist = bb.getMax()[0] - q[0];
  Vector<double, 2> direction {1.0, 0.0};

  // Uncomment this to get a reasonably fast method
  if(bb.getMax()[1] - q[1] < edge_dist)
  {
    direction = Vector<T, 2> {0.0, 1.0};
    edge_dist = bb.getMax()[1] - q[1];
  }
  if(q[0] - bb.getMin()[0] < edge_dist)
  {
    direction = Vector<T, 2> {-1.0, 0.0};
    edge_dist = q[0] - bb.getMin()[0];
  }
  if(q[1] - bb.getMin()[1] < edge_dist)
  {
    direction = Vector<T, 2> {0.0, -1.0};
  }

  // Get a ray extending in a random direction
  Ray<T, 2> ray(q, direction);
  for(auto& c : cs)
  {
    inter += detail::ray_casting_bezier_clipping(c, ray, nevals, edge_tol);
  }

  return inter;
}

/*!
 * \brief Computes the generalized winding number for a 2D curved polygon
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The CurvedPolygon object
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
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

template <typename T>
double winding_number_approxogon_memoized(
  Point<T, 2> q,
  const axom::Array<std::vector<std::vector<detail::BezierCurveMemo<T>>>>& array_memos,
  axom::Array<
    axom::FlatMap<std::pair<int, int>, detail::BezierCurveMemo<T>, detail::PairHash>>&
    hash_memos,
  Polygon<T, 2>& temp_approxogon,
  std::stack<std::pair<int, int>>& curve_stack,
  double edge_tol = 1e-8,
  double linear_tol = 1e-8)
{
  double wn = 0;

  for(int ci = 0; ci < array_memos.size(); ++ci)
  {
    auto& the_curve = array_memos[ci][0][0].curve;

    if(!the_curve.boundingBox().contains(q))
    {
      wn += detail::linear_winding_number(q,
                                          the_curve[0],
                                          the_curve[the_curve.getOrder()],
                                          edge_tol);
    }
    else
    {
      //temp_approxogon.addVertex(the_curve[0]);
      detail::winding_number_adaptive_linear_memoized(q,
                                                      array_memos[ci],
                                                      hash_memos[ci],
                                                      edge_tol,
                                                      linear_tol,
                                                      temp_approxogon,
                                                      curve_stack,
                                                      wn);

      temp_approxogon.addVertex(the_curve[the_curve.getOrder()]);

      wn += detail::approxogon_winding_number(q, temp_approxogon, edge_tol);

      //desmos_print(c);
      //desmos_print(q);
      //desmos_print(temp_approxogon);

      temp_approxogon.clear();
      //int xx = 0;
    }
  }

  return wn;
}

template <typename T>
double winding_number_approxogon(const Point<T, 2>& q,
                                 const axom::Array<BezierCurve<T, 2>>& carray,
                                 Polygon<T, 2>& temp_approxogon,
                                 int& num_evals,
                                 int& max_depth,
                                 double edge_tol = 1e-8,
                                 double linear_tol = 1e-8)
{
  double wn = 0;

  for(int i = 0; i < carray.size(); i++)
  {
    // Check exterior bounding box
    if(!carray[i].boundingBox().contains(q))
    {
      wn += detail::linear_winding_number(q,
                                          carray[i][0],
                                          carray[i][carray[i].getOrder()],
                                          edge_tol);
    }
    else
    {
      //wn += detail::curve_winding_number_recursive(q,
      //                                             carray[i],
      //                                             false,
      //                                             dummy_val,
      //                                             edge_tol);
      temp_approxogon.addVertex(carray[i][0]);
      detail::winding_number_adaptive_linear(q,
                                             carray[i],
                                             false,
                                             num_evals,
                                             edge_tol,
                                             linear_tol,
                                             temp_approxogon,
                                             wn);
      temp_approxogon.addVertex(carray[i][carray[i].getOrder()]);
      max_depth =
        axom::utilities::max(max_depth, temp_approxogon.numVertices() - 2);
      wn += detail::approxogon_winding_number(q, temp_approxogon, edge_tol);
      temp_approxogon.clear();
    }
  }

  return wn;
}

template <typename T>
double winding_number(const Point<T, 2>& q,
                      const axom::Array<BezierCurve<T, 2>>& carray,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  double ret_val = 0.0;
  int dummy_val = 0;
  for(int i = 0; i < carray.size(); i++)
  {
    ret_val += detail::curve_winding_number_recursive(q,
                                                      carray[i],
                                                      false,
                                                      dummy_val,
                                                      edge_tol);
  }

  return ret_val;
}
//@}

//@{
//! @name Winding number operations between 3D points and primitives

/*!
 * \brief Computes the solid angle winding number for a 3D triangle
 *
 * \param [in] query The query point to test
 * \param [in] tri The 3D Triangle object
 * \param [in] isOnFace An optional return parameter if the point is on the triangle
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
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

  const double denom = a_norm * b_norm * c_norm  //
    + a_norm * b.dot(c) + b_norm * a.dot(c) + c_norm * a.dot(b);

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
  bool isOnFace = false;
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
  bool isOnFace = false;
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

#ifdef AXOM_USE_MFEM

/*
 * \brief Computes the solid angle winding number for a Bezier patch
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
 * \return double The generalized winding number.
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
