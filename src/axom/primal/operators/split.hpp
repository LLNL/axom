// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file split.hpp
 *
 * \brief Consists of a set of methods to split a primal primitive into
 *        (a collection of) another primal primitive
 */

#ifndef AXOM_PRIMAL_SPLIT_HPP_
#define AXOM_PRIMAL_SPLIT_HPP_

#include "axom/core/Array.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/is_convex.hpp"

#include "axom/slic.hpp"

// perhaps #include split_impl.hpp here

namespace axom
{
namespace primal
{
/*!
 * \brief Splits an Octahedron into eight Tetrahedrons
 *
 * \tparam Tp the coordinate type, such double or float
 * \tparam NDIMS the number of spatial dimensions (must be 3).
 * \param [in] oct The Octahedron to split
 * \param [out] out The \a axom::Array of Tetrahedron objects; the fragments of
 *              oct are appended to out.
 *
 * \pre NDIMS == 3
 *
 * The tets are produced by putting a vertex at the centroid of the oct
 * and drawing an edge from each vertex to the centroid.
 * 
 */
template <typename Tp, int NDIMS = 3>
void split(const Octahedron<Tp, NDIMS>& oct,
           axom::Array<Tetrahedron<Tp, NDIMS>>& out)
{
  // Implemented for three dimensions
  SLIC_ASSERT(NDIMS == 3);

  // Type aliases
  using NumArray = NumericArray<Tp, NDIMS>;
  using Oct = Octahedron<Tp, NDIMS>;
  using Tet = Tetrahedron<Tp, NDIMS>;

  // Step 1: Find the centroid
  NumArray c;  // ctor fills with 0
  for(int i = 0; i < Oct::NUM_OCT_VERTS; ++i)
  {
    c += oct[i].array();
  }
  c = c / (double)Oct::NUM_OCT_VERTS;
  typename Oct::PointType C(c);

  // Step 2: Now store the new tets.  The documentation for the Octahedron class
  // shows how the points are arranged to imply the faces.

  // clang-format off
  enum OctVerts {P, Q, R, S, T, U};
  // clang-format on

  out.push_back(Tet(oct[P], oct[R], oct[Q], C));
  out.push_back(Tet(oct[Q], oct[R], oct[S], C));
  out.push_back(Tet(oct[R], oct[T], oct[S], C));
  out.push_back(Tet(oct[S], oct[T], oct[U], C));
  out.push_back(Tet(oct[T], oct[P], oct[U], C));
  out.push_back(Tet(oct[U], oct[P], oct[Q], C));
  out.push_back(Tet(oct[P], oct[T], oct[R], C));
  out.push_back(Tet(oct[Q], oct[S], oct[U], C));
};

/*!
 * \brief Splits a BezierCurve object into convex Bezier Curves
 *
 * \param [in] c BezierCurve to split
 * \param [out] out The \a axom::Array of BezierCurve objects; a set of convex
 *              components of c are appended to out.
 *
 * \pre NDIMS == 2
 *
 * Uses a bisection method, splitting the curve recursively until each section
 * is convex.
 * 
 */
template <typename Tp, int NDIMS = 2>
void split_to_convex(const BezierCurve<Tp, NDIMS>& c,
                     axom::Array<BezierCurve<Tp, NDIMS>>& out)
{
  // Implemented for two dimensions
  SLIC_ASSERT(NDIMS == 2);

  using Poly = Polygon<Tp, NDIMS>;
  using Bez = BezierCurve<Tp, NDIMS>;

  if(is_convex(Poly(c.getControlPoints())))
    out.push_back(Bez(c));
  else
  {
    Bez c1, c2;
    c.split(0.5, c1, c2);
    split_to_convex(c1, out);
    split_to_convex(c2, out);
  }
}

template <typename Tp, int NDIMS = 2>
void split_to_convex(const CurvedPolygon<Tp, NDIMS>& cpoly,
                     axom::Array<BezierCurve<Tp, NDIMS>>& out)
{
  // Implemented for two dimensions
  SLIC_ASSERT(NDIMS == 2);

  using Poly = Polygon<Tp, NDIMS>;
  using Bez = BezierCurve<Tp, NDIMS>;

  for(int i = 0; i < cpoly.numEdges(); i++) split_to_convex(cpoly[i], out);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_SPLIT_HPP_
