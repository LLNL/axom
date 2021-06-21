// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_IMPL_HPP_
#define AXOM_SPIN_BVH_IMPL_HPP_

#include "axom/core/Types.hpp"              // fixed bitwidth types
#include "axom/core/execution/for_all.hpp"  // for generic for_all()
#include "axom/core/memory_management.hpp"  // for memory functions
#include "axom/core/numerics/floating_point_limits.hpp"  // floating_point_limits
#include "axom/core/utilities/AnnotationMacros.hpp"      // for annotations

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/primal/operators/intersect.hpp"

// slic includes
#include "axom/slic/interface/slic.hpp"  // for SLIC macros

// linear bvh includes
#include "axom/spin/internal/linear_bvh/build_radix_tree.hpp"
#include "axom/spin/internal/linear_bvh/bvh_traverse.hpp"
#include "axom/spin/internal/linear_bvh/bvh_vtkio.hpp"
#include "axom/spin/internal/linear_bvh/BVHData.hpp"
#include "axom/spin/internal/linear_bvh/emit_bvh.hpp"
#include "axom/spin/internal/linear_bvh/QueryAccessor.hpp"

// RAJA includes
#include "RAJA/RAJA.hpp"

// C/C++ includes
#include <fstream>  // for std::ofstream
#include <sstream>  // for std::ostringstream
#include <string>   // for std::string
#include <cstring>  // for memcpy

namespace axom
{
namespace spin
{
template <typename FloatType>
using floating_point_limits = axom::numerics::floating_point_limits<FloatType>;

//------------------------------------------------------------------------------
//  BVH IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  PRIVATE HELPER METHOD IMPLEMENTATION
//------------------------------------------------------------------------------
namespace lbvh = internal::linear_bvh;

/*!
 * \def BVH_PREDICATE
 *
 * \brief Macro that defines a traversal predicate functor.
 *
 * \param _predicateName the name of the predicate, e.g., `leftPredicate`
 * \param _p the primitive type
 * \param _bb the AABB for the current BVH bin
 *
 * \note The BVH_PREDICATE may be instantiated within or outside a kernel.
 *
 * \note This macro is intended to be used internally by the BVH implementation.
 */
#define BVH_PREDICATE(_predicateName, _p, _bb) \
  auto _predicateName = [=] AXOM_HOST_DEVICE(_p, _bb) -> bool

/*!
 * \def BVH_LEAF_ACTION
 *
 * \brief Macro that defines the leaf action functor for a BVH traversal.
 *
 * \param _funcName the nem of the functor, .e.g, `leafAction`
 * \param _node the BVH node ID
 * \param _leafNodes the leafNodes
 *
 * \note The BVH_LEAF_ACTION macro must be called within a kernel.
 *
 * \note This macro is intended to be used internally by the BVH implementation.
 */
#define BVH_LEAF_ACTION(_funcName, _node, _leafNodes) \
  auto _funcName = [&](_node, _leafNodes) -> void

namespace
{
/*!
 * \brief Performs a traversal to count the candidates for each query point.
 *
 * \param [in] binCheck traversal predicate functor for bin check.
 * \param [in] inner_nodes array of vec4s for the BVH inner nodes.
 * \param [in] leaf_nodes array of BVH leaf node indices
 * \param [in] N the number of user-supplied query points
 * \param [out] counts array of candidate counts for each query point.
 * \param [in] x user-supplied array of x-coordinates
 * \param [in] y user-supplied array of y-coordinates
 * \param [in] z user-supplied array of z-coordinates
 *
 * \return total_count the total count of candidates for all query points.
 */
template <int NDIMS, typename ExecSpace, typename Predicate, typename FloatType>
IndexType bvh_get_counts(Predicate&& binCheck,
                         const primal::BoundingBox<FloatType, NDIMS>* inner_nodes,
                         const int32* inner_node_children,
                         const int32* leaf_nodes,
                         IndexType N,
                         IndexType* counts,
                         const FloatType* x,
                         const FloatType* y,
                         const FloatType* z) noexcept
{
  AXOM_PERF_MARK_FUNCTION("bvh_get_pointcounts");

  // sanity checks
  SLIC_ASSERT(inner_nodes != nullptr);
  SLIC_ASSERT(inner_node_children != nullptr);
  SLIC_ASSERT(leaf_nodes != nullptr);
  SLIC_ERROR_IF(counts == nullptr, "supplied null pointer for counts!");
  SLIC_ERROR_IF(x == nullptr, "supplied null pointer for x-coordinates!");
  SLIC_ERROR_IF(y == nullptr, "supplied null pointer for y-coordinates!");
  SLIC_ERROR_IF((z == nullptr && NDIMS == 3),
                "supplied null pointer for z-coordinates!");

  // STEP 1: count number of candidates for each query point
  using reduce_pol = typename axom::execution_space<ExecSpace>::reduce_policy;
  RAJA::ReduceSum<reduce_pol, IndexType> total_count(0);

  using QueryAccessor = lbvh::QueryAccessor<NDIMS, FloatType>;
  for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(IndexType i) {
      int32 count = 0;
      primal::Point<FloatType, NDIMS> point;
      QueryAccessor::getPoint(point, i, x, y, z);

      BVH_LEAF_ACTION(leafAction,
                      int32 AXOM_NOT_USED(current_node),
                      const int32* AXOM_NOT_USED(leaf_nodes))
      {
        count++;
      };

      lbvh::bvh_traverse(inner_nodes,
                         inner_node_children,
                         leaf_nodes,
                         point,
                         binCheck,
                         leafAction);

      counts[i] = count;
      total_count += count;
    });

  return (total_count.get());
}

/*!
 * \brief Performs a traversal to count the number of candidates for each ray.
 *
 * \param [in] binCheck traversal predicate functor for bin check.
 * \param [in] inner_nodes array of vec4s for the BVH inner nodes.
 * \param [in] leaf_nodes array of BVH leaf node indices
 * \param [in] N the number of user-supplied rays in query.
 * \param [out] counts array of length N with candidate counts for each ray.
 * \param [in] x0 array of length N with ray source point x-coordinates.
 * \param [in] nx array of length N with ray normal x-components.
 * \param [in] y0 array of length N with ray source point y-coordinates.
 * \param [in] ny array of length N with ray normal y-components.
 * \param [in] z0 array of length N with ray source point z-coordinates.
 * \param [in] nz array of length N with ray normal z-components.
 *
 * \return total_count the aggregate number of candidates for all rays.
 */
template <int NDIMS, typename ExecSpace, typename Predicate, typename FloatType>
IndexType bvh_get_raycounts(Predicate&& binCheck,
                            const primal::BoundingBox<FloatType, NDIMS>* inner_nodes,
                            const int32* inner_node_children,
                            const int32* leaf_nodes,
                            IndexType N,
                            IndexType* counts,
                            const FloatType* x0,
                            const FloatType* nx,
                            const FloatType* y0,
                            const FloatType* ny,
                            const FloatType* z0,
                            const FloatType* nz) noexcept
{
  AXOM_PERF_MARK_FUNCTION("bvh_get_raycounts");

  // sanity checks
  SLIC_ASSERT(inner_nodes != nullptr);
  SLIC_ASSERT(inner_node_children != nullptr);
  SLIC_ASSERT(leaf_nodes != nullptr);
  SLIC_ERROR_IF(counts == nullptr, "supplied null pointer for counts!");
  SLIC_ERROR_IF(x0 == nullptr, "ray source x-coordinates array is null!");
  SLIC_ERROR_IF(nx == nullptr, "ray normal x-components array is null!");
  SLIC_ERROR_IF(y0 == nullptr, "ray source y-coordinates array is null!");
  SLIC_ERROR_IF(ny == nullptr, "ray normal y-components array is null!");
  SLIC_ERROR_IF((z0 == nullptr && NDIMS == 3),
                "ray source z-coordinates array is null!");
  SLIC_ERROR_IF((nz == nullptr && NDIMS == 3),
                "ray normal z-components is null!");

  // STEP 1: count number of candidates for each query point
  using reduce_pol = typename axom::execution_space<ExecSpace>::reduce_policy;
  RAJA::ReduceSum<reduce_pol, IndexType> total_count(0);

  using QueryAccessor = lbvh::QueryAccessor<NDIMS, FloatType>;
  for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(IndexType i) {
      int32 count = 0;

      primal::Point<FloatType, NDIMS> origin;
      primal::Vector<FloatType, NDIMS> direction;
      QueryAccessor::getPoint(origin, i, x0, y0, z0);
      QueryAccessor::getPoint(direction, i, nx, ny, nz);
      primal::Ray<FloatType, NDIMS> ray {origin, direction};

      BVH_LEAF_ACTION(leafAction,
                      int32 AXOM_NOT_USED(current_node),
                      const int32* AXOM_NOT_USED(leaf_nodes))
      {
        count++;
      };

      lbvh::bvh_traverse<NDIMS>(inner_nodes,
                                inner_node_children,
                                leaf_nodes,
                                ray,
                                binCheck,
                                leafAction);

      counts[i] = count;
      total_count += count;
    });

  return (total_count.get());
}

/*!
 * \brief Performs a traversal to count the number of candidates for each
 *  bounding box.
 *
 * \param [in] binCheck traversal predicate functor for bin check.
 * \param [in] inner_nodes array of vec4s for the BVH inner nodes.
 * \param [in] leaf_nodes array of BVH leaf node indices
 * \param [in] N the number of user-supplied bounding boxes in query.
 * \param [out] counts array of length N with candidate counts for each box.
 * \param [in] xmin array of x-coordinate of the lower bounding box corner
 * \param [in] xmax array of x-coordinate of the upper bounding box corner
 * \param [in] ymin array of y-coordinate of the lower bounding box corner
 * \param [in] ymax array of y-coordinate of the upper bounding box corner
 * \param [in] zmin array of z-coordinate of the lower bounding box corner
 * \param [in] zmax array of z-coordinate of the upper bounding box corner
 *
 * \return total_count the aggregate number of candidates for all bounding
 *  boxes.
 */
template <int NDIMS, typename ExecSpace, typename Predicate, typename FloatType>
IndexType bvh_get_boxcounts(Predicate&& binCheck,
                            const primal::BoundingBox<FloatType, NDIMS>* inner_nodes,
                            const int32* inner_node_children,
                            const int32* leaf_nodes,
                            IndexType N,
                            IndexType* counts,
                            const FloatType* xmin,
                            const FloatType* xmax,
                            const FloatType* ymin,
                            const FloatType* ymax,
                            const FloatType* zmin,
                            const FloatType* zmax) noexcept
{
  AXOM_PERF_MARK_FUNCTION("bvh_get_boxcounts");

  // sanity checks
  SLIC_ASSERT(inner_nodes != nullptr);
  SLIC_ASSERT(inner_node_children != nullptr);
  SLIC_ASSERT(leaf_nodes != nullptr);
  SLIC_ERROR_IF(counts == nullptr, "supplied null pointer for counts!");
  SLIC_ERROR_IF(xmin == nullptr,
                "bounding box lower x-coordinates array is null!");
  SLIC_ERROR_IF(xmax == nullptr,
                "bounding box upper x-coordinates array is null!");
  SLIC_ERROR_IF(ymin == nullptr,
                "bounding box lower y-coordinates array is null!");
  SLIC_ERROR_IF(ymax == nullptr,
                "bounding box upper y-coordinates array is null!");
  SLIC_ERROR_IF((zmin == nullptr && NDIMS == 3),
                "bounding box lower z-coordinates array is null!");
  SLIC_ERROR_IF((zmax == nullptr && NDIMS == 3),
                "bounding box upper z-coordinates array is null!");

  // STEP 1: count number of candidates for each query point
  using reduce_pol = typename axom::execution_space<ExecSpace>::reduce_policy;
  RAJA::ReduceSum<reduce_pol, IndexType> total_count(0);

  using QueryAccessor = lbvh::QueryAccessor<NDIMS, FloatType>;
  for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(IndexType i) {
      int32 count = 0;
      primal::BoundingBox<FloatType, NDIMS> box;
      QueryAccessor::getBoundingBox(box, i, xmin, xmax, ymin, ymax, zmin, zmax);

      BVH_LEAF_ACTION(leafAction,
                      int32 AXOM_NOT_USED(current_node),
                      const int32* AXOM_NOT_USED(leaf_nodes))
      {
        count++;
      };

      lbvh::bvh_traverse(inner_nodes,
                         inner_node_children,
                         leaf_nodes,
                         box,
                         binCheck,
                         leafAction);

      counts[i] = count;
      total_count += count;
    });

  return (total_count.get());
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  PUBLIC API IMPLEMENTATION
//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType>
BVH<NDIMS, ExecSpace, FloatType>::BVH(const FloatType* boxes,
                                      IndexType numItems,
                                      int allocatorID)
  : m_AllocatorID(allocatorID)
  , m_Tolernace(floating_point_limits<FloatType>::epsilon())
  , m_scaleFactor(DEFAULT_SCALE_FACTOR)
  , m_numItems(numItems)
  , m_boxes(boxes)
{ }

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType>
BVH<NDIMS, ExecSpace, FloatType>::~BVH()
{
  m_bvh.deallocate();
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType>
int BVH<NDIMS, ExecSpace, FloatType>::build()
{
  AXOM_PERF_MARK_FUNCTION("BVH::build");

  // STEP 1: Handle case when user supplied a single bounding box
  int numBoxes = m_numItems;
  FloatType* boxesptr = nullptr;
  if(m_numItems == 1)
  {
    numBoxes = 2;
    constexpr int32 M = NDIMS * 2;  // number of entries for one box
    const int N = numBoxes * M;     // number of entries for N boxes
    boxesptr = axom::allocate<FloatType>(N, m_AllocatorID);

    const FloatType* myboxes = m_boxes;

    // copy first box and add a fake 2nd box
    for_all<ExecSpace>(
      N,
      AXOM_LAMBDA(IndexType i) { boxesptr[i] = (i < M) ? myboxes[i] : 0.0; });

  }  // END if single item
  else
  {
    boxesptr = const_cast<FloatType*>(m_boxes);
  }

  // STEP 2: Build a RadixTree consisting of the bounding boxes, sorted
  // by their corresponding morton code.
  lbvh::RadixTree<FloatType, NDIMS> radix_tree;
  primal::BoundingBox<FloatType, NDIMS> global_bounds;
  lbvh::build_radix_tree<ExecSpace>(boxesptr,
                                    numBoxes,
                                    global_bounds,
                                    radix_tree,
                                    m_scaleFactor,
                                    m_AllocatorID);

  // STEP 3: emit the BVH data-structure from the radix tree
  m_bvh.m_bounds = global_bounds;
  m_bvh.allocate(numBoxes, m_AllocatorID);

  // STEP 4: emit the BVH
  lbvh::emit_bvh<ExecSpace>(radix_tree, m_bvh);

  radix_tree.deallocate();

  // STEP 5: deallocate boxesptr if user supplied a single box
  if(m_numItems == 1)
  {
    SLIC_ASSERT(boxesptr != m_boxes);
    axom::deallocate(boxesptr);
  }

  return BVH_BUILD_OK;
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType>
void BVH<NDIMS, ExecSpace, FloatType>::getBounds(FloatType* min,
                                                 FloatType* max) const
{
  SLIC_ASSERT(min != nullptr);
  SLIC_ASSERT(max != nullptr);
  for(int idim = 0; idim < NDIMS; idim++)
  {
    min[idim] = m_bvh.m_bounds.getMin()[idim];
    max[idim] = m_bvh.m_bounds.getMax()[idim];
  }
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType>
void BVH<NDIMS, ExecSpace, FloatType>::findPoints(IndexType* offsets,
                                                  IndexType* counts,
                                                  IndexType*& candidates,
                                                  IndexType numPts,
                                                  const FloatType* x,
                                                  const FloatType* y,
                                                  const FloatType* z) const
{
  AXOM_PERF_MARK_FUNCTION("BVH::findPoints");

  SLIC_ASSERT(offsets != nullptr);
  SLIC_ASSERT(counts != nullptr);
  SLIC_ASSERT(candidates == nullptr);
  SLIC_ASSERT(x != nullptr);
  SLIC_ASSERT(y != nullptr);

  using PointType = primal::Point<FloatType, NDIMS>;
  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;
  using QueryAccessor = lbvh::QueryAccessor<NDIMS, FloatType>;

  // STEP 1: Grab BVH pointers
  const BoundingBoxType* inner_nodes = m_bvh.m_inner_nodes;
  const int32* inner_node_children = m_bvh.m_inner_node_children;
  const int32* leaf_nodes = m_bvh.m_leaf_nodes;
  SLIC_ASSERT(inner_nodes != nullptr);
  SLIC_ASSERT(inner_node_children != nullptr);
  SLIC_ASSERT(leaf_nodes != nullptr);

  // STEP 2: define traversal predicates
  BVH_PREDICATE(predicate, const PointType& p, const BoundingBoxType& bb)
  {
    return bb.contains(p);
  };

  // STEP 3: get counts
  int total_count = 0;
  AXOM_PERF_MARK_SECTION(
    "PASS[1]:count_traversal",
    total_count = bvh_get_counts<NDIMS, ExecSpace>(predicate,
                                                   inner_nodes,
                                                   inner_node_children,
                                                   leaf_nodes,
                                                   numPts,
                                                   counts,
                                                   x,
                                                   y,
                                                   z););

  using exec_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  AXOM_PERF_MARK_SECTION(
    "exclusive_scan",
    RAJA::exclusive_scan<exec_policy>(counts,
                                      counts + numPts,
                                      offsets,
                                      RAJA::operators::plus<IndexType> {}););

  AXOM_PERF_MARK_SECTION(
    "allocate_candidates",
    IndexType total_candidates = static_cast<IndexType>(total_count);
    candidates = axom::allocate<IndexType>(total_candidates, m_AllocatorID););

  // STEP 4: fill in candidates for each point
  AXOM_PERF_MARK_SECTION(
    "PASS[2]:fill_traversal",
    for_all<ExecSpace>(
      numPts,
      AXOM_LAMBDA(IndexType i) {
        int32 offset = offsets[i];

        PointType point;
        QueryAccessor::getPoint(point, i, x, y, z);

        BVH_LEAF_ACTION(leafAction, int32 current_node, const int32* leaf_nodes)
        {
          candidates[offset] = leaf_nodes[current_node];
          offset++;
        };

        lbvh::bvh_traverse(inner_nodes,
                           inner_node_children,
                           leaf_nodes,
                           point,
                           predicate,
                           leafAction);
      }););
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType>
void BVH<NDIMS, ExecSpace, FloatType>::findRays(IndexType* offsets,
                                                IndexType* counts,
                                                IndexType*& candidates,
                                                IndexType numRays,
                                                const FloatType* x0,
                                                const FloatType* nx,
                                                const FloatType* y0,
                                                const FloatType* ny,
                                                const FloatType* z0,
                                                const FloatType* nz) const
{
  AXOM_PERF_MARK_FUNCTION("BVH::findRays");

  SLIC_ASSERT(offsets != nullptr);
  SLIC_ASSERT(counts != nullptr);
  SLIC_ASSERT(candidates == nullptr);
  SLIC_ASSERT(x0 != nullptr);
  SLIC_ASSERT(nx != nullptr);
  SLIC_ASSERT(y0 != nullptr);
  SLIC_ASSERT(ny != nullptr);

  const FloatType TOL = m_Tolernace;

  using RayType = primal::Ray<FloatType, NDIMS>;
  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;
  using QueryAccessor = lbvh::QueryAccessor<NDIMS, FloatType>;

  // STEP 1: Grab BVH pointers
  const BoundingBoxType* inner_nodes = m_bvh.m_inner_nodes;
  const int32* inner_node_children = m_bvh.m_inner_node_children;
  const int32* leaf_nodes = m_bvh.m_leaf_nodes;
  SLIC_ASSERT(inner_nodes != nullptr);
  SLIC_ASSERT(inner_node_children != nullptr);
  SLIC_ASSERT(leaf_nodes != nullptr);

  // STEP 2: define traversal predicates
  BVH_PREDICATE(predicate, const RayType& r, const BoundingBoxType& bb)
  {
    primal::Point<FloatType, NDIMS> tmp;
    return primal::detail::intersect_ray(r, bb, tmp, TOL);
  };

  // STEP 3: get counts
  int total_count = 0;
  AXOM_PERF_MARK_SECTION(
    "PASS[1]:count_traversal",
    total_count = bvh_get_raycounts<NDIMS, ExecSpace>(predicate,
                                                      inner_nodes,
                                                      inner_node_children,
                                                      leaf_nodes,
                                                      numRays,
                                                      counts,
                                                      x0,
                                                      nx,
                                                      y0,
                                                      ny,
                                                      z0,
                                                      nz););

  using exec_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  AXOM_PERF_MARK_SECTION(
    "exclusive_scan",
    RAJA::exclusive_scan<exec_policy>(counts,
                                      counts + numRays,
                                      offsets,
                                      RAJA::operators::plus<IndexType> {}););

  AXOM_PERF_MARK_SECTION(
    "allocate_candidates",
    IndexType total_candidates = static_cast<IndexType>(total_count);
    candidates = axom::allocate<IndexType>(total_candidates, m_AllocatorID););

  // STEP 4: fill in candidates for each point
  AXOM_PERF_MARK_SECTION(
    "PASS[2}:fill_traversal",
    for_all<ExecSpace>(
      numRays,
      AXOM_LAMBDA(IndexType i) {
        int32 offset = offsets[i];

        typename RayType::PointType origin;
        typename RayType::VectorType direction;
        QueryAccessor::getPoint(origin, i, x0, y0, z0);
        QueryAccessor::getPoint(direction, i, nx, ny, nz);
        RayType ray {origin, direction};

        BVH_LEAF_ACTION(leafAction, int32 current_node, const int32* leaf_nodes)
        {
          candidates[offset] = leaf_nodes[current_node];
          offset++;
        };

        lbvh::bvh_traverse(inner_nodes,
                           inner_node_children,
                           leaf_nodes,
                           ray,
                           predicate,
                           leafAction);
      }););
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType>
void BVH<NDIMS, ExecSpace, FloatType>::findBoundingBoxes(IndexType* offsets,
                                                         IndexType* counts,
                                                         IndexType*& candidates,
                                                         IndexType numBoxes,
                                                         const FloatType* xmin,
                                                         const FloatType* xmax,
                                                         const FloatType* ymin,
                                                         const FloatType* ymax,
                                                         const FloatType* zmin,
                                                         const FloatType* zmax) const
{
  AXOM_PERF_MARK_FUNCTION("BVH::findBoundingBoxes");

  SLIC_ASSERT(offsets != nullptr);
  SLIC_ASSERT(counts != nullptr);
  SLIC_ASSERT(candidates == nullptr);
  SLIC_ASSERT(xmin != nullptr);
  SLIC_ASSERT(xmax != nullptr);
  SLIC_ASSERT(ymin != nullptr);
  SLIC_ASSERT(ymax != nullptr);

  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;
  using QueryAccessor = lbvh::QueryAccessor<NDIMS, FloatType>;

  // STEP 1: Grab BVH pointers
  const BoundingBoxType* inner_nodes = m_bvh.m_inner_nodes;
  const int32* inner_node_children = m_bvh.m_inner_node_children;
  const int32* leaf_nodes = m_bvh.m_leaf_nodes;
  SLIC_ASSERT(inner_nodes != nullptr);
  SLIC_ASSERT(inner_node_children != nullptr);
  SLIC_ASSERT(leaf_nodes != nullptr);

  // STEP 2: define traversal predicates
  BVH_PREDICATE(predicate, const BoundingBoxType& bb1, const BoundingBoxType& bb2)
  {
    return bb1.intersectsWith(bb2);
  };

  // STEP 3: get counts
  int total_count = 0;
  AXOM_PERF_MARK_SECTION(
    "PASS[1]:count_traversal",
    total_count = bvh_get_boxcounts<NDIMS, ExecSpace>(predicate,
                                                      inner_nodes,
                                                      inner_node_children,
                                                      leaf_nodes,
                                                      numBoxes,
                                                      counts,
                                                      xmin,
                                                      xmax,
                                                      ymin,
                                                      ymax,
                                                      zmin,
                                                      zmax););

  using exec_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  AXOM_PERF_MARK_SECTION(
    "exclusive_scan",
    RAJA::exclusive_scan<exec_policy>(counts,
                                      counts + numBoxes,
                                      offsets,
                                      RAJA::operators::plus<IndexType> {}););

  AXOM_PERF_MARK_SECTION(
    "allocate_candidates",
    IndexType total_candidates = static_cast<IndexType>(total_count);
    candidates = axom::allocate<IndexType>(total_candidates, m_AllocatorID););

  // STEP 4: fill in candidates for each bounding box
  AXOM_PERF_MARK_SECTION(
    "PASS[2}:fill_traversal",
    for_all<ExecSpace>(
      numBoxes,
      AXOM_LAMBDA(IndexType i) {
        int32 offset = offsets[i];

        BoundingBoxType box;
        QueryAccessor::getBoundingBox(box, i, xmin, xmax, ymin, ymax, zmin, zmax);

        BVH_LEAF_ACTION(leafAction, int32 current_node, const int32* leaf_nodes)
        {
          candidates[offset] = leaf_nodes[current_node];
          offset++;
        };

        lbvh::bvh_traverse(inner_nodes,
                           inner_node_children,
                           leaf_nodes,
                           box,
                           predicate,
                           leafAction);
      }););
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace, typename FloatType>
void BVH<NDIMS, ExecSpace, FloatType>::writeVtkFile(const std::string& fileName) const
{
  std::ostringstream nodes;
  std::ostringstream cells;
  std::ostringstream levels;

  // STEP 0: Write VTK header
  std::ofstream ofs;
  ofs.open(fileName.c_str());
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " BVHTree \n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  // STEP 1: write root
  int32 numPoints = 0;
  int32 numBins = 0;
  lbvh::write_root(m_bvh.m_bounds, numPoints, numBins, nodes, cells, levels);

  // STEP 2: traverse the BVH and dump each bin
  constexpr int32 ROOT = 0;
  lbvh::write_recursive<FloatType, NDIMS>(m_bvh.m_inner_nodes,
                                          m_bvh.m_inner_node_children,
                                          ROOT,
                                          1,
                                          numPoints,
                                          numBins,
                                          nodes,
                                          cells,
                                          levels);

  // STEP 3: write nodes
  ofs << "POINTS " << numPoints << " double\n";
  ofs << nodes.str() << std::endl;

  // STEP 4: write cells
  const int32 nnodes = (NDIMS == 2) ? 4 : 8;
  ofs << "CELLS " << numBins << " " << numBins * (nnodes + 1) << std::endl;
  ofs << cells.str() << std::endl;

  // STEP 5: write cell types
  ofs << "CELL_TYPES " << numBins << std::endl;
  const int32 cellType = (NDIMS == 2) ? 9 : 12;
  for(int32 i = 0; i < numBins; ++i)
  {
    ofs << cellType << std::endl;
  }

  // STEP 6: dump level information
  ofs << "CELL_DATA " << numBins << std::endl;
  ofs << "SCALARS level int\n";
  ofs << "LOOKUP_TABLE default\n";
  ofs << levels.str() << std::endl;
  ofs << std::endl;

  // STEP 7: close file
  ofs.close();
}

#undef BVH_PREDICATE
#undef BVH_LEAF_ACTION

} /* namespace spin */
} /* namespace axom */

#endif /* AXOM_SPIN_BVH_IMPL_HPP_ */
