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

// Algortihm adapted from [Maa 99]
template <typename T>
bool my_isConvex(const Polygon<T, 2>& poly)
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
      if(res1 == orientation(poly[0], seg)) return false;
    }
    else
    {
      if(res1 == orientation(poly[n], seg)) return false;
    }
  }

  return true;
}

// Check if point is interior to polygon
template <typename T>
bool my_containsPoint(const Polygon<T, 2>& poly,
                      const Point<T, 2>& p,
                      const double EPS = 1e-8)
{
  Ray<T, 2> the_ray(p, Vector<T, 2>({1, 0}));
  double ray_param = -1, seg_param = -1;
  int n = poly.numVertices() - 1;

  Segment<T, 2> the_seg(poly[n], poly[0]);

  int num_intersects = intersect_ray(the_ray, the_seg, ray_param, seg_param, EPS);
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
    if(this_intersect == true && ray_param < EPS)
      return true;
    else
      ++num_intersects;
  }

  // Return true if num_intersects is odd
  return (num_intersects % 2) == 1;
}

// Close the bezier curve, compute winding number of that closure.
template <typename T>
double closure_winding_number(BezierCurve<T, 2> c, Point2D q)
{
  Triangle<double, 2> T(q, c[0], c[c.getOrder()]);
  return 0.5 * M_1_PI * T.angle(0) * ((T.signedArea() < 0) ? -1.0 : 1.0);
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

// Recursive winding number function
template <typename T>
double adaptive_winding_number(const BezierCurve<T, 2>& c,
                               const Point2D& q,
                               int qnodes,
                               double itol,
                               double ltol,
                               int depth = 0)
{
  //auto winding_func = get_winding_func(q);
  double cl_winding_num = closure_winding_number(c, q);

  Polygon<T, 2> controlPolygon = c.getControlPolygon();

  if(my_isConvex(controlPolygon) && !my_containsPoint(controlPolygon, q, ltol))
  {
    //std::cout << "LOL using advanced algorithm" << std::endl;
    return -cl_winding_num;
  }

  /*
  double winding_num = evaluate_line_integral(c, winding_func, qnodes);
  double closed_winding = winding_num + cl_winding_num;

  // If the value is close to an integer, we know our quadrature is accurate.
  //  Return the value to machine precision.
  if(std::abs(std::round(closed_winding) - closed_winding) < itol)
  {
    return std::round( closed_winding ) - cl_winding_num;
  }
  // Or, if the curve is linear and our query point is close,
  //  then we only need to know which side of the curve we are on.
  //  Serves as the base case for our recursion
  */
  else if(c.isLinear(ltol))
  {
    /*
    std::cout << std::endl
              << "needed linear tolerance at " << depth << " depth crimging rn "
              << std::endl;
    */
    return -cl_winding_num;
  }
  // Otherwise, our quadrature didn't give us a good enough answer, so we try again
  else
  {
    BezierCurve<double, 2> c1, c2;
    c.split(0.5, c1, c2);

    // clang-format off
    double wn1 = adaptive_winding_number(c1, q, qnodes, itol, ltol, depth + 1);
    double wn2 = adaptive_winding_number(c2, q, qnodes, itol, ltol, depth + 1);
    
    /*
    std::cout << "----" << depth << "----" << std::endl;
    std::cout << wn1 << ": ";
    for( int i = 0; i < c1.getOrder() + 1; i++ )
      std::cout << axom::fmt::format("({0:.10f}, {1:.10f}), ", c1[i][0], c1[i][1]);

    std::cout << std::endl << wn2 << ": ";
    for( int i = 0; i < c2.getOrder() + 1; i++ )
      std::cout << axom::fmt::format("({0:.10f}, {1:.10f}), ", c2[i][0], c2[i][1]);
    std::cout << std::endl;
    // clang-format on
    */

    return wn1 + wn2;
  }
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif