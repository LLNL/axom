// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_SCATTERED_INTERPOLATION_H_
#define QUEST_SCATTERED_INTERPOLATION_H_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"

#include "axom/fmt.hpp"

#include "conduit_blueprint.hpp"

#include <list>
#include <vector>
#include <set>
#include <cstdlib>
#include <cmath>

namespace
{
/// Helper function to extract the dimension from the coordinate values group
/// of a mesh blueprint coordset
inline int extractDimension(const conduit::Node& values_node)
{
  SLIC_ASSERT(values_node.has_child("x"));
  return values_node.has_child("z") ? 3 : (values_node.has_child("y") ? 2 : 1);
}

/// Helper function to extract the number of points from the coordinate values group
/// of a mesh blueprint coordset
inline int extractSize(const conduit::Node& values_node)
{
  SLIC_ASSERT(values_node.has_child("x"));
  return values_node["x"].dtype().number_of_elements();
}

/**
 * \brief Utility function to create an axom::ArrayView over the array
 * of native types stored by a conduit::Node
 */
template <typename T>
inline axom::ArrayView<T> ArrayView_from_Node(conduit::Node& node, int sz)
{
  T* ptr = node.value();
  return axom::ArrayView<T>(ptr, sz);
}

/**
 * \brief Template specialization of ArrayView_from_Node for Point<double,2>
 *
 * \warning Assumes the underlying data is an MCArray with stride 2 access
 */
template <>
inline axom::ArrayView<axom::primal::Point<double, 2>> ArrayView_from_Node(
  conduit::Node& node,
  int sz)
{
  using PointType = axom::primal::Point<double, 2>;

  PointType* ptr = static_cast<PointType*>(node.data_ptr());
  return axom::ArrayView<PointType>(ptr, sz);
}

/**
 * \brief Template specialization of ArrayView_from_Node for Point<double,3>
 *
 * \warning Assumes the underlying data is an MCArray with stride 3 access
 */
template <>
inline axom::ArrayView<axom::primal::Point<double, 3>> ArrayView_from_Node(
  conduit::Node& node,
  int sz)
{
  using PointType = axom::primal::Point<double, 3>;

  PointType* ptr = static_cast<PointType*>(node.data_ptr());
  return axom::ArrayView<PointType>(ptr, sz);
}

}  // namespace

namespace axom
{
namespace quest
{
/**
 * \brief A class to perform scattered data interpolation at arbitrary points
 * over an input point set
 *
 * The class uses linear interpolation over a Delaunay triangulation of the point set.
 */
template <int NDIMS = 3>
class ScatteredInterpolation
{
public:
  static constexpr int DIM = NDIMS;
  using DelaunayTriangulation = Delaunay<DIM>;
  using PointType = typename DelaunayTriangulation::PointType;
  using BoundingBoxType = typename DelaunayTriangulation::BoundingBox;

  void buildTriangulation(conduit::Node& mesh_node, const std::string& coordset)
  {
    // Perform some simple error checking
    SLIC_ASSERT(this->isValidBlueprint(mesh_node));

    auto valuesPath = fmt::format("coordsets/{}/values", coordset);
    SLIC_ASSERT(mesh_node.has_path(valuesPath));
    auto& values = mesh_node[valuesPath];

    SLIC_ASSERT(::extractDimension(values) == DIM);
    const int npts = ::extractSize(values);

    // Assume initially that coords are interleaved
    auto coords = ::ArrayView_from_Node<PointType>(values["x"], npts);

    // Compute the bounding box
    BoundingBoxType bb;
    for(const auto& pt : coords)
    {
      bb.addPoint(pt);
    }
    bb.scale(1.5);

    m_delaunay.initializeBoundary(bb);
    for(const auto& pt : coords)
    {
      m_delaunay.insertPoint(pt);
    }

    m_delaunay.removeBoundary();

    m_delaunay.writeToVTKFile(fmt::format("delaunay_{}d.vtk", DIM));
  }

private:
  /// Check validity of blueprint group
  bool isValidBlueprint(const conduit::Node& mesh_node) const
  {
    bool success = true;
    conduit::Node info;
    if(!conduit::blueprint::verify("mesh", mesh_node, info))
    {
      SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
      success = false;
    }

    return success;
  }

private:
  DelaunayTriangulation m_delaunay;
};

template <int NDIMS>
constexpr int ScatteredInterpolation<NDIMS>::DIM;

}  // namespace quest
}  // namespace axom

#endif  // QUEST_SCATTERED_INTERPOLATION_H_