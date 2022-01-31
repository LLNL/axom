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
  void setObjectMesh(const sidre::Group* meshGroup)
  {
    // Perform some simple error checking
    SLIC_ASSERT(meshGroup != nullptr);
    SLIC_ASSERT(meshGroup->hasGroup("coordsets/coords/values"));
    SLIC_ASSERT(meshGroup->hasView("coordsets/coords/values/x"));

    // Extract the dimension and number of points
    const auto* valuesGroup = meshGroup->getGroup("coordsets/coords/values");
    const int dim =
      valuesGroup->hasView("z") ? 3 : (valuesGroup->hasChildView("y") ? 2 : 1);
    const int N = valuesGroup->getView("x")->getNumElements();
    SLIC_ASSERT(dim == NDIMS);

    // Extract pointers to the coordinate data
    // Note: The following assumes that the coordinate arrays have stride 1
    // TODO: Generalize to support other strides
    primal::ZipIndexable<PointType> zipPts {
      {valuesGroup->hasView("x")
         ? static_cast<const double*>(valuesGroup->getView("x")->getVoidPtr())
         : nullptr,
       valuesGroup->hasView("y")
         ? static_cast<const double*>(valuesGroup->getView("y")->getVoidPtr())
         : nullptr,
       valuesGroup->hasView("z")
         ? static_cast<const double*>(valuesGroup->getView("z")->getVoidPtr())
         : nullptr}};

    // Copy the data into an array of primal points
    PointArray pts(0, N);
    for(int i = 0; i < N; ++i)
    {
      pts.push_back(zipPts[i]);
    }

    // Optionally, output some diagnostic information about the input points
    if(m_isVerbose)
    {
      SLIC_INFO("Points on object:");
      for(auto i : slam::PositionSet<int>(pts.size()))
      {
        const double mag = sqrt(pts[i][0] * pts[i][0] + pts[i][1] * pts[i][1]);
        SLIC_INFO(axom::fmt::format("\t{}: {} -- {}", i, pts[i], mag));
      }
    }

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

  void computeClosestPoints(const PointArray& queryPts,
                            axom::Array<axom::IndexType>& cpIndexes) const
  {
    SLIC_ASSERT(!queryPts.empty());

    const int nPts = queryPts.size();

    /// Create an ArrayView in ExecSpace that is compatible with cpIndexes
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

    cpIndexes = query_inds;
  }

  const PointArray& points() const { return m_points; }

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
