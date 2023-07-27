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
#include "axom/primal/geometry/BezierPatch.hpp"

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
  for(int i = 0; i < Oct::NUM_VERTS; ++i)
  {
    c += oct[i].array();
  }
  c = c / static_cast<double>(Oct::NUM_VERTS);
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
 * \brief Splits a BezierPatch object into SuperConvex subpatches
 *
 * \param [in] bPatch BezierPatch to split
 * \param [out] out The \a axom::Array of BezierPatch objects; a set of SuperConvex
 *              components of c are appended to out.
 *
 *
 * Uses a bisection method, splitting the patch recursively until each section
 * is SuperConvex (convex + shallow).
 * 
 */
template <typename Tp>
void split_to_valid(const BezierPatch<Tp>& bPatch,
                    axom::Array<BezierPatch<Tp>>& out)
{
  using Poly = Polygon<Tp, NDIMS>;
  using Patch = BezierPatch<Tp>;

  const int ord_u = bPatch.getOrder_u();
  const int ord_v = bPatch.getOrder_v();
  Poly control_slice();

  bool is_valid = true;

  // Check the slices that are fixed in u
  for(int p = 0; is_valid && p <= ord_u; ++p)
  {
    control_slice.clear();
    for(int q = 0; q <= ord_v; ++q)
    {
      control_slice.addVertex(bPatch(p, q));
    }

    is_valid = is_convex(control_slice) && is_shallow(control_slice);
  }

  // Check the slices that are fixed in v
  for(int q = 0; is_valid && q <= ord_v; ++q)
  {
    control_slice.clear();
    for(int p = 0; p <= ord_u; ++p)
    {
      control_slice.addVertex(bPatch(p, q));
    }
    is_valid = is_convex(control_slice) && is_shallow(control_slice);
  }

  if(is_valid)
  {
    out.push_back(Patch(bPatch));
  }
  else
  {
    Patch p1, p2, p3, p4;
    bPatch.split(0.5, p1, p2, p3, p4);
    split_to_valid(p1, out);
    split_to_valid(p2, out);
    split_to_valid(p3, out);
    split_to_valid(p4, out);
  }
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_SPLIT_HPP_
