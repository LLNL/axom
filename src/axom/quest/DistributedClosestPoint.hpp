// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DISTRIBUTED_CLOSEST_POINT_H_
#define QUEST_DISTRIBUTED_CLOSEST_POINT_H_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/sidre.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"

#include "axom/fmt.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mpi.hpp"

#include <list>
#include <vector>
#include <set>
#include <cstdlib>
#include <cmath>

namespace axom
{
namespace quest
{
template <int NDIMS = 2, typename ExecSpace = axom::SEQ_EXEC>
class DistributedClosestPoint
{
public:
  // TODO: generalize to 3D
  static_assert(NDIMS == 2,
                "DistributedClosestPoint only currently supports 2D");

  static constexpr int DIM = NDIMS;
  using PointType = primal::Point<double, DIM>;
  using BoxType = primal::BoundingBox<double, DIM>;
  using PointArray = axom::Array<PointType>;
  using BoxArray = axom::Array<BoxType>;
  using BVHTreeType = spin::BVH<DIM, ExecSpace>;

private:
  struct MinCandidate
  {
    /// Squared distance to query point
    double minSqDist {numerics::floating_point_limits<double>::max()};
    /// Index within mesh of closest element
    int minElem;
  };

public:
  DistributedClosestPoint(
    int allocatorID = axom::execution_space<ExecSpace>::allocatorID())
    : m_allocatorID(allocatorID)
  { }

public:  // Query properties
  void setVerbosity(bool isVerbose) { m_isVerbose = isVerbose; }

public:
  /** 
   * Utility function to set the array of points
   *
   * \param [in] meshGroup The root group of a mesh blueprint
   * \note This function currently supports mesh blueprints with the "point" topology
   */
  void setObjectMesh(sidre::Group* meshGroup, const std::string& coordset)
  {
    // Perform some simple error checking
    SLIC_ASSERT(this->isValidBlueprint(meshGroup));

    // Extract the dimension and number of points from the coordinate values group
    auto valuesPath = axom::fmt::format("coordsets/{}/values", coordset);
    SLIC_ASSERT(meshGroup->hasGroup(valuesPath));
    auto* valuesGroup = meshGroup->getGroup(valuesPath);
    const int dim = extractDimension(valuesGroup);
    const int N = extractSize(valuesGroup);
    SLIC_ASSERT(dim == NDIMS);

    // Extract pointers to the coordinate data
    // Note: The following assumes the coordinates are contiguous with stride NDIMS
    // TODO: Generalize to support other strides

    // Copy the data into the point array of primal points
    PointArray pts(N, N);
    const std::size_t nbytes = sizeof(double) * dim * N;
    axom::copy(pts.data(),
               valuesGroup->getView("x")->getBuffer()->getVoidPtr(),
               nbytes);

    m_points = PointArray(pts, m_allocatorID);  // copy point array to ExecSpace
  }

  bool generateBVHTree()
  {
    const int npts = m_points.size();
    axom::Array<BoxType> boxesArray(npts, npts, m_allocatorID);
    auto boxesView = boxesArray.view();
    axom::for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(axom::IndexType i) { boxesView[i] = BoxType {m_points[i]}; });

    // Build bounding volume hierarchy
    m_bvh.setAllocatorID(m_allocatorID);
    int result = m_bvh.initialize(boxesView, npts);

    return (result == spin::BVH_BUILD_OK);
  }

  void computeClosestPoints(sidre::Group* meshGroup,
                            const std::string& coordset) const
  {
    // Perform some simple error checking
    SLIC_ASSERT(this->isValidBlueprint(meshGroup));

    // Extract the dimension and number of points from the coordinate values group
    auto valuesPath = axom::fmt::format("coordsets/{}/values", coordset);
    SLIC_ASSERT(meshGroup->hasGroup(valuesPath));
    auto* valuesGroup = meshGroup->getGroup(valuesPath);
    const int dim = extractDimension(valuesGroup);
    const int nPts = extractSize(valuesGroup);
    SLIC_ASSERT(dim == NDIMS);

    // Extract the points array
    auto* data = static_cast<PointType*>(valuesGroup->getView("x")->getVoidPtr());
    axom::ArrayView<PointType> queryPts(data, nPts);

    // Extract the cp_index array
    std::string cp_index_path = "fields/cp_index/values";
    SLIC_ASSERT_MSG(
      meshGroup->hasView(cp_index_path),
      "Input to `computeClosestPoint()` must have a field names `cp_index`");
    sidre::Array<axom::IndexType> cpIndexes(meshGroup->getView(cp_index_path));

    // Create an ArrayView in ExecSpace that is compatible with cpIndexes
    // TODO: Avoid copying (here and at the end) if both are on the host
    axom::Array<axom::IndexType> cp_idx(nPts, nPts, m_allocatorID);
    auto query_inds = cp_idx.view();

    /// Create an ArrayView in ExecSpace that is compatible with queryPts
    PointArray execPoints(nPts, nPts, m_allocatorID);
    execPoints = queryPts;
    auto query_pts = execPoints.view();

    // Get a device-useable iterator
    auto it = m_bvh.getTraverser();

    using axom::primal::squared_distance;
    using int32 = axom::int32;

    AXOM_PERF_MARK_SECTION(
      "ComputeClosestPoints",
      axom::for_all<ExecSpace>(
        nPts,
        AXOM_LAMBDA(int32 idx) mutable {
          PointType qpt = query_pts[idx];

          MinCandidate curr_min {};

          auto searchMinDist = [&](int32 current_node, const int32* leaf_nodes) {
            const int candidate_idx = leaf_nodes[current_node];
            const PointType candidate_pt = m_points[candidate_idx];
            const double sq_dist = squared_distance(qpt, candidate_pt);

            if(sq_dist < curr_min.minSqDist)
            {
              curr_min.minSqDist = sq_dist;
              curr_min.minElem = candidate_idx;
            }
          };

          auto traversePredicate = [&](const PointType& p,
                                       const BoxType& bb) -> bool {
            return squared_distance(p, bb) <= curr_min.minSqDist;
          };

          // Traverse the tree, searching for the point with minimum distance.
          it.traverse_tree(qpt, searchMinDist, traversePredicate);

          query_inds[idx] = curr_min.minElem;
        }););

    axom::copy(cpIndexes.data(),
               query_inds.data(),
               cpIndexes.size() * sizeof(axom::IndexType));
  }

  const PointArray& points() const { return m_points; }

private:
  // Check validity of blueprint group
  bool isValidBlueprint(const sidre::Group* bpGroup) const
  {
    if(bpGroup == nullptr)
    {
      return false;
    }

    conduit::Node mesh_node;
    bpGroup->createNativeLayout(mesh_node);

    bool success = true;
    conduit::Node info;
    if(!conduit::blueprint::mpi::verify("mesh", mesh_node, info, MPI_COMM_WORLD))
    {
      SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
      success = false;
    }
    // else
    // {
    //   SLIC_INFO("Valid blueprint for particle mesh: \n" << info.to_yaml());
    // }

    return success;
  }

  /// Helper function to extract the dimension from the coordinate values group
  int extractDimension(sidre::Group* valuesGroup) const
  {
    SLIC_ASSERT(valuesGroup->hasView("x"));
    return valuesGroup->hasView("z") ? 3 : (valuesGroup->hasView("y") ? 2 : 1);
  }

  /// Helper function to extract the number of points from the coordinate values group
  int extractSize(sidre::Group* valuesGroup) const
  {
    SLIC_ASSERT(valuesGroup->hasView("x"));
    return valuesGroup->getView("x")->getNumElements();
  }

private:
  PointArray m_points;
  BoxArray m_boxes;
  BVHTreeType m_bvh;

  int m_allocatorID;
  bool m_isVerbose {false};
};

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DISTRIBUTED_CLOSEST_POINT_H_
