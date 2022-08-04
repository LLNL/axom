// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file is_linear.hpp
 *
 * \brief Consists of a method to evaluate if a BezierCurve object is approximately linear.
 *
 * Does this by ensuring that the control nodes are within tol of a segment connecting
 * the endpoints of the curve. This ensures that the entire curve is also within
 * tol of this segment.
 * 
 */

#ifndef PRIMAL_IS_LINEAR_HPP_
#define PRIMAL_IS_LINEAR_HPP_

// Axom includes
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"

namespace axom
{
namespace primal
{

template <typename T, int NDIMS>
bool is_linear(BezierCurve<T, NDIMS> curve, double tol = 1E-8)
{
  const int ord = curve.getOrder();
  if(ord <= 1)
  {
    return true;
  }

  Segment<T, NDIMS> seg(curve[0], curve[ord]);
  double sqDist = 0.0;
  for(int p = 1; p < ord && sqDist < tol; ++p)  // check interior control points
  {
    sqDist += squared_distance(curve[p], seg);
  }
  return (sqDist < tol);
}

}  // end namespace primal
}  // end namespace axom

#endif