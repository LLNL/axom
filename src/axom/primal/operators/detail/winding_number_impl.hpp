// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_WINDING_NUMBER_IMPL_HPP_
#define PRIMAL_WINDING_NUMBER_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/in_polygon.hpp"
#include "axom/primal/operators/is_convex.hpp"
#include "axom/primal/operators/squared_distance.hpp"
#include "axom/primal/operators/intersect.hpp"

// C++ includes
#include <math.h>
#include <unordered_map>
#include <utility>
#include <stack>

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

namespace axom
{
namespace primal
{
namespace detail
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
 * Colinear points return 0 for their winding number.
 *
 * \return double The winding number
 */
template <typename T>
double linear_winding_number(const Point<T, 2>& q,
                             const Point<T, 2>& c0,
                             const Point<T, 2>& c1,
                             double edge_tol)
{
  Vector<T, 2> V1(q, c0);
  Vector<T, 2> V2(q, c1);

  // clang-format off
  // Measures the signed area of the triangle with vertices q, c0, c1
  double tri_area = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                                V1[1] - V2[1], V2[1]);
  // clang-format on

  // Compute distance from line connecting endpoints to query
  if(tri_area * tri_area <= edge_tol * edge_tol * (V1 - V2).squared_norm())
  {
    return 0;
  }

  // Compute signed angle between vectors
  double dotprod = axom::utilities::clampVal(
    Vector<T, 2>::dot_product(V1.unitVector(), V2.unitVector()),
    -1.0,
    1.0);

  return 0.5 * M_1_PI * acos(dotprod) * ((tri_area > 0) ? 1 : -1);
}

/*!
 * \brief Directly compute the winding number at either endpoint of a 
 *        Bezier curve with a convex control polygon
 *
 * \param [in] q The query point
 * \param [in] c The BezierCurve object to compute the winding number along
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance for isNearlyZero
 * \pre Control polygon for c must be convex
 * \pre The query point must be on one of the endpoints
 *
 * The winding number for a Bezier curve with a convex control polygon is
 * given by the signed angle between the tangent vector at that endpoint and
 * the vector in the direction of the other endpoint. 
 * 
 * The query can be located on both endpoints if it is closed, in which case
 * the angle is that between the tangent lines at both endpoints
 * 
 * \return double The winding number
 */
template <typename T>
double convex_endpoint_winding_number(const Point<T, 2>& q,
                                      const BezierCurve<T, 2>& c,
                                      double edge_tol)
{
  const int ord = c.getOrder();
  if(ord == 1)
  {
    return 0;
  }

  double edge_tol_sq = edge_tol * edge_tol;

  // Verify that the shape is convex, and that the query point is at an endpoint
  SLIC_ASSERT(is_convex(Polygon<T, 2>(c.getControlPoints()), PRIMAL_TINY));
  SLIC_ASSERT((squared_distance(q, c[0]) <= edge_tol_sq) ||
              (squared_distance(q, c[ord]) <= edge_tol_sq));

  int idx;

  // Need to find vectors that subtend the entire curve.
  //   We must ignore duplicate nodes
  for(idx = 0; idx <= ord; ++idx)
  {
    if(squared_distance(q, c[idx]) > edge_tol_sq)
    {
      break;
    }
  }
  Vector<T, 2> V1(q, c[idx]);

  for(idx = ord; idx >= 0; --idx)
  {
    if(squared_distance(q, c[idx]) > edge_tol_sq)
    {
      break;
    }
  }
  Vector<T, 2> V2(q, c[idx]);

  // clang-format off
  // Measures the signed area of the triangle spanned by V1 and V2
  double tri_area = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                                V1[1] - V2[1], V2[1]);
  // clang-format on

  // This means the bounding vectors are anti-parallel.
  //  Parallel tangents can't happen with nontrivial convex control polygons
  if((ord > 3) && axom::utilities::isNearlyEqual(tri_area, 0.0, PRIMAL_TINY))
  {
    for(int i = 1; i < ord; ++i)
    {
      // Need to find the first non-parallel control node
      V2 = Vector<T, 2>(q, c[i]);

      // clang-format off
      tri_area = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                             V1[1] - V2[1], V2[1]);
      // clang-format on

      // Because we are convex, a single non-collinear vertex tells us the orientation
      if(!axom::utilities::isNearlyEqual(tri_area, 0.0, PRIMAL_TINY))
      {
        return (tri_area > 0) ? 0.5 : -0.5;
      }
    }

    // If all vectors are parallel, the curve is linear and return 0
    return 0;
  }

  // Compute signed angle between vectors
  double dotprod = axom::utilities::clampVal(
    Vector<T, 2>::dot_product(V1.unitVector(), V2.unitVector()),
    -1.0,
    1.0);
  return 0.5 * M_1_PI * acos(dotprod) * ((tri_area > 0) ? 1 : -1);
}

/*!
 * \brief Recursively compute the winding number for a query point with respect
 *        to a single Bezier curve.
 *
 * \param [in] q The query point at which to compute winding number
 * \param [in] c The BezierCurve object along which to compute the winding number
 * \param [in] isConvexControlPolygon Boolean flag if the input Bezier curve 
                                      is already convex
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance for nonphysical distances, used in
 *   isLinear, isNearlyZero, in_polygon, is_convex
 * 
 * Use a recursive algorithm that checks if the query point is exterior to
 * some convex shape containing the Bezier curve, in which case we have a direct 
 * formula for the winding number. If not, we bisect our curve and run the algorithm on 
 * each half. Use the proximity of the query point to endpoints and approximate
 * linearity of the Bezier curve as base cases.
 * 
 * \return double The winding number.
 */
template <typename T>
double curve_winding_number_recursive(const Point<T, 2>& q,
                                      const BezierCurve<T, 2>& c,
                                      bool isConvexControlPolygon,
                                      int& nevals,
                                      double edge_tol = 1e-8)
{
  const int ord = c.getOrder();
  if(ord <= 0)
  {
    return 0.0;  // Catch degenerate cases
  }

  // If q is outside a convex shape that contains the entire curve, the winding
  //   number for the shape connected at the endpoints with straight lines is zero.
  //   We then subtract the contribution of this line segment.

  // Simplest convex shape containing c is its bounding box
  BoundingBox<T, 2> bBox(c.boundingBox());
  if(!bBox.contains(q))
  {
    return 0.0 - linear_winding_number(q, c[ord], c[0], edge_tol);
  }

  // Use linearity as base case for recursion.
  if(c.isLinear(edge_tol))
  {
    return linear_winding_number(q, c[0], c[ord], edge_tol);
  }

  // Check if our control polygon is convex.
  //  If so, all subsequent control polygons will be convex as well
  Polygon<T, 2> controlPolygon(c.getControlPoints());
  const bool includeBoundary = true;
  const bool useNonzeroRule = true;

  if(!isConvexControlPolygon)
  {
    isConvexControlPolygon = is_convex(controlPolygon, PRIMAL_TINY);
  }
  else  // Formulas for winding number only work if shape is convex
  {
    // Bezier curves are always contained in their convex control polygon
    if(!in_polygon(q, controlPolygon, includeBoundary, useNonzeroRule, PRIMAL_TINY))
    {
      return 0.0 - linear_winding_number(q, c[ord], c[0], edge_tol);
    }

    // If the query point is at either endpoint, use direct formula
    if((squared_distance(q, c[0]) <= edge_tol * edge_tol) ||
       (squared_distance(q, c[ord]) <= edge_tol * edge_tol))
    {
      return convex_endpoint_winding_number(q, c, edge_tol);
    }
  }

  // Recursively split curve until query is outside some known convex region
  BezierCurve<T, 2> c1, c2;
  c.split(0.5, c1, c2);
  nevals += 1;

  return curve_winding_number_recursive(q,
                                        c1,
                                        isConvexControlPolygon,
                                        nevals,
                                        edge_tol) +
    curve_winding_number_recursive(q, c2, isConvexControlPolygon, nevals, edge_tol);
}

template <typename T>
double approxogon_winding_number(const Point<T, 2>& q,
                                 const Polygon<T, 2>& approxogon,
                                 double edge_tol)
{
  bool isOnEdge = false;
  const int n = approxogon.numVertices();

  // Catch some early edge cases
  if(n <= 1) return 0.0;

  int integer_wn = winding_number(q, approxogon, isOnEdge, false, PRIMAL_TINY);
  double closure_wn =
    linear_winding_number(q, approxogon[n - 1], approxogon[0], edge_tol);

  if(isOnEdge)
  {
    // If on edge, can't use integer winding number. Have to compute it
    //  through winding number of each edge

    //std::cout << "on edge" << std::endl;

    // Main leg is oriented in the reverse
    double wn = -closure_wn;
    for(int i = 1; i < n; ++i)
    {
      wn += detail::linear_winding_number(q,
                                          approxogon[i - 1],
                                          approxogon[i],
                                          edge_tol);
    }

    return wn;
  }

  return integer_wn - closure_wn;
}

template <typename T>
struct BezierCurveMemo
{
  bool isConvexControlPolygon;
  BezierCurve<T, 2> curve;
};

struct PairHash
{
  using result_type = std::size_t;

  size_t operator()(const std::pair<int, int>& p) const
  {
    // Combine hashes of the two integers
    size_t hash1 = std::hash<int>()(p.first);
    size_t hash2 = std::hash<int>()(p.second);
    return hash1 ^ (hash2 << 1);  // Simple combining function
  }
};

//bool operator==(const std::pair<int, int>& lhs, const std::pair<int, int>& rhs)
//{
//  return lhs.first == rhs.first && lhs.second == rhs.second;
//}

template <typename T>
void winding_number_adaptive_linear_memoized(
  Point<T, 2> q,
  const std::vector<std::vector<BezierCurveMemo<T>>>& array_memo,
  axom::FlatMap<std::pair<int, int>, BezierCurveMemo<T>, PairHash>& hash_memo,
  double edge_tol,
  double linear_tol,
  Polygon<T, 2>& approxogon,
  std::stack<std::pair<int, int>>& curve_stack,
  double& end_wn)
{
  constexpr int LEVEL = 3;
  const int ord = array_memo[0][0].curve.getOrder();

  //std::stack<std::pair<int, int>> curve_stack;
  curve_stack.push(std::make_pair(0, 0));

  while(!curve_stack.empty())
  {
    std::pair<int, int> the_pair = curve_stack.top();
    curve_stack.pop();

    //std::cout << the_pair.first << std::endl;

    //std::cout << the_pair.first << " " << the_pair.second << std::endl;
    BezierCurveMemo<T> curve_memo;
    if(the_pair.first < LEVEL)
      curve_memo = array_memo[the_pair.first][the_pair.second];
    else
      curve_memo = hash_memo[the_pair];

    //if(curve_memo.isLinear)
    //{
    //  approxogon.addVertex(curve_memo.curve[0]);
    //  continue;
    //}

    BoundingBox<T, 2> bBox(curve_memo.curve.boundingBox());
    if(!bBox.contains(q))
    {
      approxogon.addVertex(curve_memo.curve[0]);
      continue;
    }

    if(curve_memo.isConvexControlPolygon)
    {
      constexpr bool includeBoundary = true;
      constexpr bool useNonzeroRule = true;

      if(!in_polygon(q,
                     Polygon<T, 2>(curve_memo.curve.getControlPoints()),
                     includeBoundary,
                     useNonzeroRule,
                     PRIMAL_TINY))
      {
        approxogon.addVertex(curve_memo.curve[0]);

        continue;
      }

      if((squared_distance(q, curve_memo.curve[0]) <= edge_tol * edge_tol) ||
         (squared_distance(q, curve_memo.curve[ord]) <= edge_tol * edge_tol))
      {
        end_wn += approxogon_winding_number(q, approxogon, edge_tol);
        approxogon.clear();

        end_wn += convex_endpoint_winding_number(q, curve_memo.curve, edge_tol);
        approxogon.addVertex(curve_memo.curve[ord]);
        continue;
      }
    }

    auto ref1 = std::make_pair(the_pair.first + 1, 2 * the_pair.second);
    auto ref2 = std::make_pair(the_pair.first + 1, 2 * the_pair.second + 1);

    if(the_pair.first >= LEVEL - 1 && hash_memo.find(ref1) == hash_memo.end())
    {
      BezierCurve<T, 2> c1, c2;
      curve_memo.curve.split(0.5, c1, c2);

      hash_memo[ref1] = {
        //c1.isLinear(linear_tol),
        curve_memo.isConvexControlPolygon ||
          is_convex(Polygon<T, 2>(c1.getControlPoints()), PRIMAL_TINY),
        c1};

      hash_memo[ref2] = {
        //c2.isLinear(linear_tol),
        curve_memo.isConvexControlPolygon ||
          is_convex(Polygon<T, 2>(c2.getControlPoints()), PRIMAL_TINY),
        c2};
    }

    curve_stack.push(ref2);
    curve_stack.push(ref1);
  }
}

template <typename T>
void winding_number_adaptive_linear(const Point<T, 2>& q,
                                    const BezierCurve<T, 2>& c,
                                    bool isConvexControlPolygon,
                                    int& num_evals,
                                    double edge_tol,
                                    double linear_tol,
                                    Polygon<T, 2>& approxogon,
                                    double& end_wn)
{
  const int ord = c.getOrder();
  // If q is outside a convex shape that contains the entire curve, the winding
  //   number for the shape connected at the endpoints with straight lines is zero.
  //   We then subtract the contribution of this line segment.

  // Simplest convex shape containing c is its bounding box
  BoundingBox<T, 2> bBox(c.boundingBox());
  if(!bBox.contains(q))
  {
    // Don't need to bisect any further
    //desmos_print(c);
    return;
  }

  // Use linearity as base case for recursion.
  if(c.isLinear(linear_tol))
  {
    // todo: why do we need to start over if we hit this case?

    //std::cout << "is linear case" << std::endl;
    //std::cout << approxogon.numVertices() << std::endl;
    //desmos_print(c);

    // Don't need to bisect any further, but we do need
    //  to start the polygon over,
    //end_wn += approxogon_winding_number(q, approxogon, edge_tol);
    //approxogon.clear();

    // and use the direct formula for the line segment
    //end_wn += linear_winding_number(q, c[0], c[ord], edge_tol);
    //++num_evals;
    //approxogon.addVertex(c[ord]);
    return;
  }

  // Check if our control polygon is convex.
  //  If so, all subsequent control polygons will be convex as well
  Polygon<T, 2> controlPolygon(c.getControlPoints());
  const bool includeBoundary = true;
  const bool useNonzeroRule = true;

  if(!isConvexControlPolygon)
  {
    isConvexControlPolygon = is_convex(controlPolygon, PRIMAL_TINY);
  }

  // Formulas for winding number only work if shape is convex
  if(isConvexControlPolygon)
  {
    // Bezier curves are always contained in their convex control polygon
    if(!in_polygon(q, controlPolygon, includeBoundary, useNonzeroRule, PRIMAL_TINY))
    {
      // Don't need to bisect any further
      //desmos_print(c);
      return;
    }

    // If the query point is at either endpoint, use direct formula
    if(squared_distance(q, c[0]) <= edge_tol * edge_tol ||
       squared_distance(q, c[ord]) <= edge_tol * edge_tol)
    {
      //std::cout << approxogon << std::endl;

      // Need to start the polygon over,
      end_wn += approxogon_winding_number(q, approxogon, edge_tol);
      std::cout << "using endpoint formula" << std::endl;
      approxogon.clear();

      // and use the direct formula for the endpoint
      end_wn += convex_endpoint_winding_number(q, c, edge_tol);

      approxogon.addVertex(c[ord]);
      //desmos_print(c);
      return;
    }
  }

  // Do a single iteration of Newton's Method to get a better guess than 0.5
  //double split_val = 0.5;
  //Point<T, 2> eval;
  //Vector<T, 2> Dt, DtDt;

  //for(int i = 0; i < 15; ++i)
  //{
  //  c.evaluate_second_derivative(split_val, eval, Dt, DtDt);

  //  Vector<T, 2> q_eval(q, eval);
  //  if(q_eval.squared_norm() < edge_tol * edge_tol) break;

  //  split_val = axom::utilities::clampVal(
  //    split_val - Dt.dot(q_eval) / (Dt.dot(Dt) + DtDt.dot(q_eval)),
  //    0.0,
  //    1.0);
  //}

  // Recursively split curve until query is outside some known convex region
  BezierCurve<T, 2> c1, c2;
  c.split(0.5, c1, c2);

  // clang-format off
  winding_number_adaptive_linear(q, c1, isConvexControlPolygon, num_evals, edge_tol, linear_tol, approxogon, end_wn);

  //std::cout << t1 << " " << c2[0] << c.evaluate( t1 ) << std::endl;
  approxogon.addVertex(c2[0]);
  ++num_evals;
  

  winding_number_adaptive_linear(q, c2, isConvexControlPolygon, num_evals, edge_tol, linear_tol, approxogon, end_wn);
  // clang-format on

  return;
}

template <typename T>
int ray_casting_bezier_clipping(const BezierCurve<T, 2>& curve,
                                const Ray<T, 2>& ray,
                                int& nevals,
                                double tol = 1E-8)
{
  const int ord = curve.getOrder();
  const bool isRational = curve.isRational();

  // Make an array of curves to append to later
  axom::Array<BezierCurve<T, 2>> curves(0);
  curves.push_back(curve);
  int inter = 0;

  while(!curves.empty())
  {
    BezierCurve<T, 2> the_curve = curves[curves.size() - 1];
    curves.erase(curves.end() - 1);

    Point<T, 2> origin = ray.origin();
    Vector<T, 2> e1 = ray.direction().unitVector();
    Vector<T, 2> e2 = Vector<T, 2> {-e1[1], e1[0]};

    // Iterate over the control points and find the coordinates
    //  in terms of e1 and e2

    int8_t flag = ~0;
    for(int p = 0; flag && p <= ord; ++p)
    {
      auto pt_vec = Vector<T, 2>(origin, the_curve[p]);
      double alpha = pt_vec.dot(e1);
      double beta = pt_vec.dot(e2);

      if(alpha > 0 && beta > 0)  // quad 1
      {
        // Turn off bits 1, 2, 3, 5, 6
        flag &= 0x89;
      }
      else if(alpha < 0 && beta > 0)  // quad 2
      {
        // Turn off bits 0, 2, 3, 6, 7
        flag &= 0x4C;
      }
      else if(alpha < 0 && beta < 0)
      {
        // Turn off bits 0, 1, 3, 4, 7
        flag &= 0x26;
      }
      else if(alpha > 0 && beta < 0)
      {
        // Turn off bits 0, 1, 2, 4, 5
        flag &= 0x13;
      }
    }

    if(flag == 0x0)
    {
      // Can't check intersections if the curve isn't convex
      Polygon<T, 2> controlPolygon(the_curve.getControlPoints());
      if(!is_convex(controlPolygon, tol))
      {
        BezierCurve<T, 2> c1, c2;
        the_curve.split(0.5, c1, c2);
        curves.push_back(c1);
        curves.push_back(c2);
        nevals += 1;
        continue;
      }

      BoundingBox<T, 2> bb = the_curve.boundingBox();
      double size = axom::utilities::max(bb.range()[0], bb.range()[1]);

      if(size < tol)
      {
        std::cout << "is On the curve lol" << std::endl;
        return 1;
      }

      double a = e1[1];
      double b = -e1[0];
      double c = e1[0] * origin[1] - e1[1] * origin[0];
      double u, tmp;

      // Do a very sloppy way of computing intersections with the convex hull
      double umin = 1.0, umax = 0.0;
      Ray<T, 2> axis(Point<T, 2> {0.0, 0.0}, Vector<T, 2> {1.0, 0.0});

      // Find the maximum and minimum coordinates of all the intersections.
      for(double p = 0; p < ord; ++p)
      {
        double R0 = a * the_curve[p][0] + b * the_curve[p][1] + c;
        double R1 = a * the_curve[p + 1][0] + b * the_curve[p + 1][1] + c;

        if(R0 * R1 > 0) continue;

        if(the_curve.isRational())
        {
          R0 *= the_curve.getWeight(p);
          R1 *= the_curve.getWeight(p + 1);
        }

        // Compute intersection between the ray and the segment
        Segment<T, 2> seg(Point<T, 2> {p / ord, R0},
                          Point<T, 2> {(p + 1) / ord, R1});
        intersect(axis, seg, u, tmp, 0.0);

        umin = axom::utilities::min(umin, u);
        umax = axom::utilities::max(umax, u);
      }

      // Do the final check to see if the ray intersects the curve
      double R0 = a * the_curve[ord][0] + b * the_curve[ord][1] + c;
      double R1 = a * the_curve[0][0] + b * the_curve[0][1] + c;

      if(R0 * R1 < 0)
      {
        if(the_curve.isRational())
        {
          R0 *= the_curve.getWeight(ord);
          R1 *= the_curve.getWeight(0);
        }

        // Compute intersection between the ray and the segment
        Segment<T, 2> seg(Point<T, 2> {1.0, R0}, Point<T, 2> {0.0, R1});
        intersect(axis, seg, u, tmp, 0.0);

        umin = axom::utilities::min(umin, u);
        umax = axom::utilities::max(umax, u);
      }

      // Make minor numerical adjustments
      umin = 0.99 * umin;
      umax = 0.99 * umax + 0.01;

      // Heuristic to account for multiple intersections
      if(umax - umin > 0.8)
      {
        BezierCurve<T, 2> c1, c2;
        the_curve.split(0.5, c1, c2);
        curves.push_back(c1);
        curves.push_back(c2);
        nevals += 1;
      }
      else
      {
        // Now we have the min and max, so we can split the curve
        BezierCurve<T, 2> c1, c2, c3;
        the_curve.split(umin, c1, c2);
        c2.split((umax - umin) / (1.0 - umin), c2, c3);

        nevals += 2;
        curves.push_back(c1);
        curves.push_back(c2);
        curves.push_back(c3);
      }
    }
    else if(flag == 0x1)
    {
      //std::cout << "Case B: Maybe intersections" << std::endl;
      double alpha1 = Vector<T, 2>(origin, the_curve[0]).dot(e2);
      double alpha2 = Vector<T, 2>(origin, the_curve[ord]).dot(e2);

      if(alpha1 * alpha2 < 0)
      {
        (alpha1 < 0) ? inter++ : inter--;
      }
    }
    //else
    //{
    //  std::cout << "Case A: No intersections" << std::endl;
    //}
  }

  return inter;
}

/// Type to indicate orientation of singularities relative to surface
enum class SingularityAxis
{
  x,
  y,
  z,
  rotated
};

#ifdef AXOM_USE_MFEM
/*!
 * \brief Evaluates an "anti-curl" of the winding number along a curve
 *
 * \param [in] query The query point to test
 * \param [in] curve The BezierCurve object
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] npts The number of points used in each Gaussian quadrature
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * 
 * Applies a non-adaptive quadrature to a BezierCurve using one of three possible
 * "anti-curl" vector fields, the curl of each of which is equal to <x, y, z>/||x||^3.
 * With the proper "anti-curl" selected, integrating this along a closed curve is equal
 * to evaluating the winding number of a surface with that curve as the boundary.
 * 
 * \note This is only meant to be used for `winding_number<BezierPatch>()`,
 *  and the result does not make sense outside of that context.
 *
 * \return double One component of the winding number
 */
template <typename T>
double stokes_winding_number(const Point<T, 3>& query,
                             const BezierCurve<T, 3>& curve,
                             const SingularityAxis ax,
                             int npts,
                             double quad_tol)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

  double quadrature = 0.0;
  for(int q = 0; q < quad_rule.GetNPoints(); ++q)
  {
    // Get quadrature points in space (shifted by the query)
    const Vector<T, 3> node(query, curve.evaluate(quad_rule.IntPoint(q).x));
    const Vector<T, 3> node_dt(curve.dt(quad_rule.IntPoint(q).x));
    const double node_norm = node.norm();

    // Compute one of three vector field line integrals depending on
    //  the orientation of the original surface, indicated through ax.
    switch(ax)
    {
    case(SingularityAxis::x):
      quadrature += quad_rule.IntPoint(q).weight *
        (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
        (node[1] * node[1] + node[2] * node[2]) / node_norm;
      break;
    case(SingularityAxis::y):
      quadrature += quad_rule.IntPoint(q).weight *
        (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
        (node[0] * node[0] + node[2] * node[2]) / node_norm;
      break;
    case(SingularityAxis::z):
    case(SingularityAxis::rotated):
      quadrature += quad_rule.IntPoint(q).weight *
        (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
        (node[0] * node[0] + node[1] * node[1]) / node_norm;
      break;
    }
  }

  // Adaptively refine quadrature over curves if query is not far enough away
  //  from the singularity axis. If rotated, assume you need to adapt.
  bool needs_adapt = false;
  BoundingBox<T, 3> cBox(curve.boundingBox());
  Point<T, 3> centroid = cBox.getCentroid();

  switch(ax)
  {
  case(SingularityAxis::x):
    needs_adapt = (query[1] - centroid[1]) * (query[1] - centroid[1]) +
        (query[2] - centroid[2]) * (query[2] - centroid[2]) <=
      cBox.range().squared_norm();
    break;
  case(SingularityAxis::y):
    needs_adapt = (query[0] - centroid[0]) * (query[0] - centroid[0]) +
        (query[2] - centroid[2]) * (query[2] - centroid[2]) <=
      cBox.range().squared_norm();
    break;
  case(SingularityAxis::z):
    needs_adapt = (query[0] - centroid[0]) * (query[0] - centroid[0]) *
        (query[1] - centroid[1]) * (query[1] - centroid[1]) <=
      cBox.range().squared_norm();
    break;
  case(SingularityAxis::rotated):
    needs_adapt = true;
    break;
  }

  if(needs_adapt)
  {
    return stokes_winding_number_adaptive(query,
                                          curve,
                                          ax,
                                          quad_rule,
                                          quadrature,
                                          quad_tol);
  }

  return 0.25 * M_1_PI * quadrature;
}
#endif

#ifdef AXOM_USE_MFEM
/*!
 * \brief Recursively evaluates an "anti-curl" of the winding number on subcurves
 *
 * \param [in] query The query point to test
 * \param [in] curve The BezierCurve object
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] quad_rule The mfem quadrature rule object
 * \param [in] quad_coarse The integral evaluated on the original curve
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * \param [in] depth The current recursive depth
 * 
 * Recursively apply quadrature for one of three possible integrals along two halfs 
 * of a curve. The sum of this integral along the subcurves should be equal to to
 * quad_coarse. Otherwise, apply this algorithm to each half recursively.
 *
 * \note This is only meant to be used for `winding_number<BezierPatch>()`,
 *  and the result does not make sense outside of that context.
 * 
 * \return double One component of the winding number
 */
template <typename T>
double stokes_winding_number_adaptive(const Point<T, 3>& query,
                                      const BezierCurve<T, 3>& curve,
                                      const SingularityAxis ax,
                                      const mfem::IntegrationRule& quad_rule,
                                      const double quad_coarse,
                                      const double quad_tol,
                                      const int depth = 1)
{
  // Split the curve, do the quadrature over both components
  BezierCurve<T, 3> subcurves[2];
  curve.split(0.5, subcurves[0], subcurves[1]);

  double quad_fine[2] = {0.0, 0.0};
  for(int i = 0; i < 2; ++i)
  {
    for(int q = 0; q < quad_rule.GetNPoints(); ++q)
    {
      // Get quad_rulerature points in space (shifted by the query)
      const Vector<T, 3> node(query,
                              subcurves[i].evaluate(quad_rule.IntPoint(q).x));
      const Vector<T, 3> node_dt(subcurves[i].dt(quad_rule.IntPoint(q).x));
      const double node_norm = node.norm();

      // Compute one of three vector field line integrals depending on
      //  the orientation of the original surface, indicated through ax.
      switch(ax)
      {
      case(SingularityAxis::x):
        quad_fine[i] += quad_rule.IntPoint(q).weight *
          (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
          (node[1] * node[1] + node[2] * node[2]) / node_norm;
        break;
      case(SingularityAxis::y):
        quad_fine[i] += quad_rule.IntPoint(q).weight *
          (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
          (node[0] * node[0] + node[2] * node[2]) / node_norm;
        break;
      case(SingularityAxis::z):
      case(SingularityAxis::rotated):
        quad_fine[i] += quad_rule.IntPoint(q).weight *
          (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
          (node[0] * node[0] + node[1] * node[1]) / node_norm;
        break;
      }
    }
  }

  constexpr int MAX_DEPTH = 12;
  if(depth >= MAX_DEPTH ||
     axom::utilities::isNearlyEqualRelative(quad_fine[0] + quad_fine[1],
                                            quad_coarse,
                                            quad_tol,
                                            1e-10))
  {
    return 0.25 * M_1_PI * (quad_fine[0] + quad_fine[1]);
  }
  else
  {
    return stokes_winding_number_adaptive(query,
                                          subcurves[0],
                                          ax,
                                          quad_rule,
                                          quad_fine[0],
                                          quad_tol,
                                          depth + 1) +
      stokes_winding_number_adaptive(query,
                                     subcurves[1],
                                     ax,
                                     quad_rule,
                                     quad_fine[1],
                                     quad_tol,
                                     depth + 1);
  }
}
#endif

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
