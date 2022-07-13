// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file in_curved_polygon.hpp
 *
 * \brief Consists of methods that test whether a given query point is
 * inside a given curved polygon.
 *
 * Uses an adaptive winding number calculation
 */

#ifndef AXOM_PRIMAL_IN_CURVED_POLYGON_HPP_
#define AXOM_PRIMAL_IN_CURVED_POLYGON_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/primal.hpp"

#include "axom/primal/operators/detail/in_curved_polygon_impl.hpp"

// C++ includes
#include <cmath>

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#else
  #error "Primal's in/out functions for CurvedPolygon require mfem library."
#endif

namespace axom
{
namespace primal
{

/*!
 * \brief Robustly tests whether a query point lies inside curved polygon
 */
template <typename T>
inline bool in_curved_polygon(const Point<T, 2>& query,
                              const CurvedPolygon<T, 2>& cpoly,
                              double atol = 1e-5)
{
  const quadrature_nodes = 15;
  int max_depth = 0;

  double ret_val = 0.0;
  int this_depth = 0;

  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    ret_val +=
      winding_number(cpoly[i], query, quadrature_nodes, this_depth, atol);
    max_depth = std::max(max_depth, this_depth);
  }

  return !(std::round(ret_val) == 0);
}

// Base winding number function
template <typename T>
double winding_number(const CurvedPolygon<T, 2>& cpoly,
                      const Point2D& q,
                      int qnodes,
                      double int_tol = 1e-5,
                      double linear_tol = 1e-8)
{
  double ret_val = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    double this_val =
      detail::adaptive_winding_number(cpoly[i], q, qnodes, int_tol, linear_tol);
    ret_val += this_val;
  }

  return ret_val;
}

// Overload for single bezier curve
template <typename T>
double winding_number(const BezierCurve<T, 2>& c,
                      const Point2D& q,
                      int qnodes,
                      double int_tol = 1e-5,
                      double linear_tol = 1e-8)
{
  return detail::adaptive_winding_number(c, q, qnodes, int_tol, linear_tol);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_IN_CURVED_POLYGON_H_
