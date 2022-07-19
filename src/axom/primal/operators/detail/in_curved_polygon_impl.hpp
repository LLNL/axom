// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_IN_CURVED_POLYGON_IMPL_HPP_
#define PRIMAL_IN_CURVED_POLYGON_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"  // for compile-time configuration options
#include "axom/primal.hpp"

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#else
  #error "Primal's in/out functions for CurvedPolygon require mfem library."
#endif

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
namespace detail
{

// Algortihm adapted from [Maa 99].
// Determines if the polygon is convex. Will eventually
//  be placed in Polygon.hpp
template <typename T>
bool isConvex(const Polygon<T, 2>& poly)
{
  int n = poly.numVertices() - 1;
  if(n + 1 < 3) return true;  // Triangles and lines are convex

  for(int i = 1; i < n; i++)
  {
    // For each non-endpoint, check if that point and one of the endpoints
    //  are on the same side as the segment connecting the adjacent nodes
    Segment<T, 2> seg(poly[i - 1], poly[i + 1]);
    int res1 = orientation(poly[i], seg);

    // Edge case
    if(res1 == primal::ON_BOUNDARY) continue;

    if(i < n / 2)
    {
      if(res1 == orientation(poly[n], seg)) return false;
    }
    else
    {
      if(res1 == orientation(poly[0], seg)) return false;
    }
  }

  return true;
}

// Check if point is interior to polygon. Will eventually
//  be placed in Polygon.hpp
template <typename T>
bool containsPoint(const Polygon<T, 2>& poly,
                   const Point<T, 2>& p,
                   const double EPS = 1e-8)
{
  Ray<T, 2> the_ray(p, Vector<T, 2>({1, 0}));
  double ray_param = -1, seg_param = -1;
  int n = poly.numVertices() - 1;

  Segment<T, 2> the_seg(poly[n], poly[0]);

  // To avoid double counting vertices, let the segment be open at origin
  int num_intersects =
    (intersect_ray(the_ray, the_seg, ray_param, seg_param, EPS) &&
     seg_param > EPS);

  // If the ray intersects the segment at its endpoint, then the query point
  //  is on the polygon. Interpret this as "inside." Allow for tolerance
  //  consistent with intersect_ray
  if(num_intersects == 1 && ray_param < EPS) return true;

  for(int i = 0; i < n; i++)
  {
    bool this_intersect = intersect_ray(the_ray,
                                        Segment<T, 2>(poly[i], poly[i + 1]),
                                        ray_param,
                                        seg_param,
                                        EPS);

    if(this_intersect == true && seg_param > EPS)
    {
      if(ray_param < EPS)
        return true;
      else
        ++num_intersects;
    }
  }

  // Return true if num_intersects is odd
  return (num_intersects % 2) == 1;
}

// Close the bezier curve, compute winding number of that closure.
template <typename T>
double closure_winding_number(BezierCurve<T, 2> c, Point2D q)
{
  Triangle<double, 2> T(c[c.getOrder()], q, c[0]);
  return 0.5 * M_1_PI * T.angle(1) * ((T.signedArea() < 0) ? 1.0 : -1.0);
}

// Get the function to be passed into the evaluate integral function
inline std::function<Vector2D(Point2D)> get_winding_func(Point2D p)
{
  return [p](Point2D x) -> Vector2D {
    double denom =
      2 * M_PI * ((x[0] - p[0]) * (x[0] - p[0]) + (x[1] - p[1]) * (x[1] - p[1]));
    return Vector2D({-(x[1] - p[1]) / denom, (x[0] - p[0]) / denom});
  };
}

// Recursive winding number function. Assumes convex c, and bisects c
//  until q is exterior to the bounding box
template <typename T>
double adaptive_winding_number(const BezierCurve<T, 2>& c,
                               const Point2D& q,
                               int& depth,
                               double ltol)
{
  double cl_winding_num = closure_winding_number(c, q);

  // Base case for recursion
  if(c.isLinear(ltol)) return -cl_winding_num;
  
  Polygon<T, 2> controlPolygon = c.getControlPolygon();

  // Means point is far enough away to ensure we are outside the closure
  if(!containsPoint(controlPolygon, q, ltol)) return -cl_winding_num;
  // Otherwise, our quadrature didn't give us a good enough answer, so we try again
  else
  {
    depth++;

    BezierCurve<double, 2> c1, c2;
    c.split(0.5, c1, c2);

    return adaptive_winding_number(c1, q, depth, ltol) +
      adaptive_winding_number(c2, q, depth, ltol);
  }
}

// Direct winding number function. Uses approximated winding number
//   to determine which side of the curve we are on, and which
//   side of the closure we are on.
template <typename T>
double quadrature_winding_number(const BezierCurve<T, 2>& c,
                                 const Point2D& q,
                                 int qnodes)
{
  double cl_winding_num = closure_winding_number(c, q);

  if(c.getOrder() == 1)  // If the curve is linear,
    return -cl_winding_num;

  auto winding_func = get_winding_func(q);
  double winding_num_approx = evaluate_line_integral(c, winding_func, qnodes);

  // Use triangle to get orientation of control polygon.
  //  Determines the winding number of the interior
  Triangle<T, 2> tri(c[0], c[1], c[2]);
  bool bez_orientation_pos = (tri.signedArea() > 0);

  if(bez_orientation_pos)
  {
    if((winding_num_approx > 0) && (cl_winding_num > 0))
      return 1 - cl_winding_num;  // On the interior of the closure
    else
      return -cl_winding_num;
  }
  else
  {
    if((winding_num_approx < 0) && (cl_winding_num < 0))
      return -1 - cl_winding_num;  // On the interior of the closure
    else
      return -cl_winding_num;
  }
}

// Recursive winding number function that bisects the curve until
//  the quadrature + closure number is close to an integer
template <typename T>
double rounding_winding_number(const BezierCurve<T, 2>& c,
                               const Point2D& q,
                               int qnodes,
                               int& depth,
                               double itol,
                               double ltol)
{
  double cl_winding_num = closure_winding_number(c, q);

  // Base case for our recursion
  if(c.isLinear(ltol)) return -cl_winding_num;

  auto winding_func = get_winding_func(q);
  double winding_num = evaluate_line_integral(c, winding_func, qnodes);
  double closed_winding = cl_winding_num + winding_num;

  // Use triangle to get orientation of control polygon.
  //  Determines the winding number of the interior
  Triangle<T, 2> tri(c[0], c[1], c[2]);
  int interior_wn = (tri.signedArea() > 0) ? 1 : -1;

  // If the value is close to an integer, we are confident our
  //  quadrature is accurate. Return the value to machine precision.
  if(std::abs(std::round(closed_winding) - closed_winding) < itol)
  {
    if(std::round(closed_winding) == 0)
      return -cl_winding_num;
    else
      return interior_wn - cl_winding_num;
    //return std::round(closed_winding) - cl_winding_num;
  }
  // Otherwise, our quadrature didn't give us a good enough answer, so we try again
  else
  {
    BezierCurve<double, 2> c1, c2;
    c.split(0.5, c1, c2);

    return rounding_winding_number(c1, q, depth, ltol) +
      rounding_winding_number(c2, q, depth, ltol);
  }
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif