// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DISTRIBUTED_CLOSEST_POINT_H_
#define QUEST_DISTRIBUTED_CLOSEST_POINT_H_

#include "axom/config.hpp"
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

#ifndef AXOM_USE_MPI
  #error This file requires Axom to be configured with MPI
#endif
#include "mpi.h"

namespace axom
{
namespace quest
{
namespace internal
{
/** 
 * \brief Utility function to create an axom::ArrayView over the array 
 * of native types stored by a conduit::Node
 */
template <typename T>
axom::ArrayView<T> ArrayView_from_Node(conduit::Node& node, int sz)
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
axom::ArrayView<primal::Point<double, 2>> ArrayView_from_Node(conduit::Node& node,
                                                              int sz)
{
  using PointType = primal::Point<double, 2>;

  PointType* ptr = static_cast<PointType*>(node.data_ptr());
  return axom::ArrayView<PointType>(ptr, sz);
}

}  // namespace internal

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
    int minElem {-1};
    /// MPI rank of closest element
    int minRank {-1};
  };

public:
  DistributedClosestPoint(
    int allocatorID = axom::execution_space<ExecSpace>::allocatorID())
    : m_allocatorID(allocatorID)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_nranks);
  }

public:  // Query properties
  void setVerbosity(bool isVerbose) { m_isVerbose = isVerbose; }

public:
  /** 
   * Utility function to set the array of points
   *
   * \param [in] meshGroup The root group of a mesh blueprint
   * \note This function currently supports mesh blueprints with the "point" topology
   */
  void setObjectMesh(const conduit::Node& mesh_node, const std::string& coordset)
  {
    // Perform some simple error checking
    SLIC_ASSERT(this->isValidBlueprint(mesh_node));

    // Extract the dimension and number of points from the coordinate values group
    auto valuesPath = axom::fmt::format("coordsets/{}/values", coordset);
    SLIC_ASSERT(mesh_node.has_path(valuesPath));
    auto& values = mesh_node[valuesPath];

    const int dim = extractDimension(values);
    const int N = extractSize(values);
    SLIC_ASSERT(dim == NDIMS);

    // Extract pointers to the coordinate data
    // Note: The following assumes the coordinates are contiguous with stride NDIMS
    // TODO: Generalize to support other strides

    // Copy the data into the point array of primal points
    PointArray pts(N, N);
    const std::size_t nbytes = sizeof(double) * dim * N;
    axom::copy(pts.data(), values["x"].data_ptr(), nbytes);

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

  /**
   * \brief Computes the closest point within the objects for each query point
   * in the provided particle mesh, provided in the mesh blueprint rooted at \a meshGroup
   * 
   * \param meshGroup The root of a mesh blueprint for the query points
   * \param coordset The coordinate set for the query points
   * 
   * Uses the \a coordset coordinate set of the provided blueprint mesh
   * 
   * The particle mesh must contain the following fields:
   *   - cp_index: Will hold the index of the object point containing the closest point
   *   - closest_point: Will hold the position of the closest point
   * 
   * \note The current implementation assumes that the coordinates and closest_points and contiguous
   * with stride NDIMS. We intend to loosen this restriction in the future
   */
  void computeClosestPoints(conduit::Node& mesh_node,
                            const std::string& coordset) const
  {
    // Perform some simple error checking
    SLIC_ASSERT(this->isValidBlueprint(mesh_node));

    // Extract the dimension and number of points from the coordinate values group
    const auto valuesPath = axom::fmt::format("coordsets/{}/values", coordset);
    SLIC_ASSERT(mesh_node.has_path(valuesPath));
    auto& values = mesh_node[valuesPath];
    const int dim = extractDimension(values);
    const int nPts = extractSize(values);
    SLIC_ASSERT(dim == NDIMS);

    // Extract the points array
    auto queryPts = internal::ArrayView_from_Node<PointType>(values["x"], nPts);

    // Extract the cp_index array
    const std::string cp_index_path = "fields/cp_index/values";
    SLIC_ASSERT_MSG(
      mesh_node.has_path(cp_index_path),
      "Input to `computeClosestPoint()` must have a field named `cp_index`");
    auto cpIndexes =
      internal::ArrayView_from_Node<axom::IndexType>(mesh_node[cp_index_path],
                                                     nPts);

    // Extract the cp_index array
    const std::string cp_rank_path = "fields/cp_rank/values";
    SLIC_ASSERT_MSG(
      mesh_node.has_path(cp_rank_path),
      "Input to `computeClosestPoint()` must have a field named `cp_rank`");
    auto cpRanks =
      internal::ArrayView_from_Node<axom::IndexType>(mesh_node[cp_rank_path],
                                                     nPts);

    // Extract the closest points array
    const std::string cp_start_path = "fields/closest_point/values/x";
    SLIC_ASSERT_MSG(mesh_node.has_path(cp_start_path),
                    "Input to `computeClosestPoint()` must have a field named "
                    "`closest_point`");
    auto closestPts =
      internal::ArrayView_from_Node<PointType>(mesh_node[cp_start_path], nPts);

    // Create an ArrayView in ExecSpace that is compatible with cpIndexes
    // TODO: Avoid copying (here and at the end) if both are on the host
    axom::Array<axom::IndexType> cp_idx(nPts, nPts, m_allocatorID);
    auto query_inds = cp_idx.view();

    axom::Array<axom::IndexType> cp_ranks(nPts, nPts, m_allocatorID);
    auto query_ranks = cp_ranks.view();

    // Create an ArrayView in ExecSpace that is compatible with cpIndexes
    // TODO: Avoid copying (here and at the end) if both are on the host
    axom::Array<PointType> cp_pos(nPts, nPts, m_allocatorID);
    auto query_pos = cp_pos.view();

    /// Create an ArrayView in ExecSpace that is compatible with queryPts
    PointArray execPoints(nPts, nPts, m_allocatorID);
    execPoints = queryPts;
    auto query_pts = execPoints.view();

    // Get a device-useable iterator
    auto it = m_bvh.getTraverser();

    using axom::primal::squared_distance;
    using int32 = axom::int32;

    const int rank = m_rank;

    AXOM_PERF_MARK_SECTION(
      "ComputeClosestPoints",
      axom::for_all<ExecSpace>(
        nPts,
        AXOM_LAMBDA(int32 idx) mutable {
          PointType qpt = query_pts[idx];

          MinCandidate curr_min {};
          if(query_ranks[idx] >= 0)  // i.e. we've already found a candidate closest
          {
            curr_min.minSqDist = squared_distance(qpt, query_pos[idx]);
          }

          auto searchMinDist = [&](int32 current_node, const int32* leaf_nodes) {
            const int candidate_idx = leaf_nodes[current_node];
            const PointType candidate_pt = m_points[candidate_idx];
            const double sq_dist = squared_distance(qpt, candidate_pt);

            if(sq_dist < curr_min.minSqDist)
            {
              curr_min.minSqDist = sq_dist;
              curr_min.minElem = candidate_idx;
              curr_min.minRank = rank;
            }
          };

          auto traversePredicate = [&](const PointType& p,
                                       const BoxType& bb) -> bool {
            return squared_distance(p, bb) <= curr_min.minSqDist;
          };

          // Traverse the tree, searching for the point with minimum distance.
          it.traverse_tree(qpt, searchMinDist, traversePredicate);

          // If modified, update the fields that changed
          if(curr_min.minRank == rank)
          {
            query_inds[idx] = curr_min.minElem;
            query_ranks[idx] = curr_min.minRank;
            query_pos[idx] = m_points[curr_min.minElem];
          }
        }););

    axom::copy(cpIndexes.data(),
               query_inds.data(),
               cpIndexes.size() * sizeof(axom::IndexType));
    axom::copy(cpRanks.data(),
               query_ranks.data(),
               cpRanks.size() * sizeof(axom::IndexType));
    axom::copy(closestPts.data(),
               query_pos.data(),
               closestPts.size() * sizeof(PointType));

    // std::cout << fmt::format("After: closest points ({}): {}",
    //                          fmt::ptr(closestPts.data()),
    //                          closestPts)
    //           << std::endl;
  }

private:
  // Check validity of blueprint group
  bool isValidBlueprint(const conduit::Node& mesh_node) const
  {
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
  int extractDimension(const conduit::Node& values_node) const
  {
    SLIC_ASSERT(values_node.has_child("x"));
    return values_node.has_child("z") ? 3 : (values_node.has_child("y") ? 2 : 1);
  }

  /// Helper function to extract the number of points from the coordinate values group
  int extractSize(const conduit::Node& values_node) const
  {
    SLIC_ASSERT(values_node.has_child("x"));
    return values_node["x"].dtype().number_of_elements();
  }

private:
  PointArray m_points;
  BoxArray m_boxes;
  BVHTreeType m_bvh;

  int m_allocatorID;
  bool m_isVerbose {false};

  int m_rank;
  int m_nranks;
};

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DISTRIBUTED_CLOSEST_POINT_H_
