// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DISTRIBUTED_CLOSEST_POINT_H_
#define QUEST_DISTRIBUTED_CLOSEST_POINT_H_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
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
  /// Utility function to generate an array of 2D points
  void generatePoints(double radius, int numPoints)
  {
    using axom::utilities::random_real;

    // Generate in host because random_real is not yet ported to the device
    PointArray pts;
    pts.reserve(numPoints);
    for(int i = 0; i < numPoints; ++i)
    {
      const double angleInRadians = random_real(0., 2 * M_PI);
      const double rsinT = radius * std::sin(angleInRadians);
      const double rcosT = radius * std::cos(angleInRadians);

      pts.emplace_back(PointType {rcosT, rsinT});
    }

    if(m_isVerbose)
    {
      SLIC_INFO("Points on object:");
      const auto& arr = m_points;
      for(auto i : slam::PositionSet<int>(arr.size()))
      {
        const double mag = sqrt(arr[i][0] * arr[i][0] + arr[i][1] * arr[i][1]);
        SLIC_INFO(axom::fmt::format("\t{}: {} -- {}", i, arr[i], mag));
      }
    }

    m_points = PointArray(pts, m_allocatorID);  // copy to ExecSpace
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
