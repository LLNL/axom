// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/config.hpp"

#include "axom/core/StackArray.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/quest/MeshClipperStrategy.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{
namespace detail
{

using TetrahedronType = MeshClipperStrategy::TetrahedronType;
using TetVertexValues = axom::StackArray<double, TetrahedronType::NUM_VERTS>;

/*!
 * \brief Return the squared length of the longest tetrahedron edge.
 *
 * \param [in] tet the tetrahedron whose edges are measured
 * \return the maximum squared edge length
 */
AXOM_HOST_DEVICE inline double tetMaxEdgeSquared(const TetrahedronType& tet)
{
  double maxEdgeSquared = 0.0;
  for(int i = 0; i < TetrahedronType::NUM_VERTS; ++i)
  {
    for(int j = i + 1; j < TetrahedronType::NUM_VERTS; ++j)
    {
      maxEdgeSquared = axom::utilities::max(maxEdgeSquared, (tet[j] - tet[i]).squared_norm());
    }
  }
  return maxEdgeSquared;
}

/*!
 * \brief Split a positively oriented tetrahedron into eight positive children.
 *
 * The input may alias the first output tetrahedron when callers subdivide a
 * work item in place. Copying the parent vertices before writing the children
 * preserves the input in that case.
 *
 * \param [in] tet the positively oriented parent tetrahedron
 * \param [out] children the eight positively oriented child tetrahedra
 * \note \p tet may alias \p children[0].
 */
AXOM_HOST_DEVICE inline void subdivideTetByMidpoints(const TetrahedronType& tet,
                                                     TetrahedronType children[8])
{
  using Point3DType = MeshClipperStrategy::Point3DType;

  const auto v0 = tet[0];
  const auto v1 = tet[1];
  const auto v2 = tet[2];
  const auto v3 = tet[3];
  const auto m01 = Point3DType::midpoint(v0, v1);
  const auto m02 = Point3DType::midpoint(v0, v2);
  const auto m03 = Point3DType::midpoint(v0, v3);
  const auto m12 = Point3DType::midpoint(v1, v2);
  const auto m13 = Point3DType::midpoint(v1, v3);
  const auto m23 = Point3DType::midpoint(v2, v3);

  children[0] = TetrahedronType(v0, m01, m02, m03);
  children[1] = TetrahedronType(m01, v1, m12, m13);
  children[2] = TetrahedronType(m02, m12, v2, m23);
  children[3] = TetrahedronType(m03, m13, m23, v3);

  // Split the central octahedron along the m01-m23 diagonal.
  children[4] = TetrahedronType(m01, m02, m03, m23);
  children[5] = TetrahedronType(m01, m12, m02, m23);
  children[6] = TetrahedronType(m01, m13, m12, m23);
  children[7] = TetrahedronType(m01, m03, m13, m23);
}

/*!
 * \brief Compute the positive-side volume of an affine scalar field over a tetrahedron.
 *
 * The field is specified by its values at the tetrahedron vertices. Since the
 * zero isosurface is planar, its clipped-volume fraction can be computed
 * directly from the four values without constructing intersection geometry.
 *
 * \param [in] values the affine scalar-field values at the four tet vertices
 * \param [in] tetVolume the positive volume of the tetrahedron
 * \return the tetrahedron volume on the nonnegative side of the field
 * \pre \p tetVolume is nonnegative
 */
AXOM_HOST_DEVICE inline double clipTetByVertexValues(const TetVertexValues& values, double tetVolume)
{
  axom::StackArray<int, TetrahedronType::NUM_VERTS> inside;
  axom::StackArray<int, TetrahedronType::NUM_VERTS> outside;
  int insideCount = 0;
  int outsideCount = 0;
  for(int i = 0; i < TetrahedronType::NUM_VERTS; ++i)
  {
    if(values[i] >= 0.0)
    {
      inside[insideCount++] = i;
    }
    else
    {
      outside[outsideCount++] = i;
    }
  }

  if(insideCount == 0)
  {
    return 0.0;
  }
  if(insideCount == TetrahedronType::NUM_VERTS)
  {
    return tetVolume;
  }

  if(insideCount == 1 || insideCount == 3)
  {
    const bool invert = insideCount == 3;
    const int isolated = invert ? outside[0] : inside[0];
    const int* others = invert ? inside.data() : outside.data();
    double fraction = 1.0;
    for(int i = 0; i < 3; ++i)
    {
      fraction *= values[isolated] / (values[isolated] - values[others[i]]);
    }
    return invert ? tetVolume * (1.0 - fraction) : tetVolume * fraction;
  }

  const int i0 = inside[0];
  const int i1 = inside[1];
  const int o0 = outside[0];
  const int o1 = outside[1];
  const double w00 = values[i0] / (values[i0] - values[o0]);
  const double w01 = values[i0] / (values[i0] - values[o1]);
  const double w10 = values[i1] / (values[i1] - values[o0]);
  const double w11 = values[i1] / (values[i1] - values[o1]);
  const double fraction = w00 * w01 + w00 * w11 * (1.0 - w01) + w10 * w11 * (1.0 - w00);
  return tetVolume * fraction;
}

}  // namespace detail
}  // namespace experimental
}  // namespace quest
}  // namespace axom
