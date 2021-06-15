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
 * \brief Performs a traversal to find the candidates for each query primitive.
 *
 * \param [in] predicate traversal predicate functor for bin check.
 * \param [in] inner_nodes array of bounding boxes for the BVH inner nodes.
 * \param [in] inner_node_children array of children indices for the BVH inner nodes.
 * \param [in] leaf_nodes array of BVH leaf node indices
 * \param [out] offsets array of offsets into the candidate array for each query primitive
 * \param [out] counts array of candidate counts for each query primitive
 * \param [out] candidates array of the potential candidates for intersection with the BVH
 * \param [in] num_objs the number of user-supplied query primitives
 * \param [in] objs array of primitives to query against the BVH
 *
 * \return total_count the total count of candidates for all query primitives.
 */
template <int NDIMS, typename ExecSpace, typename Predicate, typename FloatType, typename PrimitiveType>
void find_candidates_impl(Predicate&& predicate,
                          const primal::BoundingBox<FloatType, NDIMS>* inner_nodes,
                          const int32* inner_node_children,
                          const int32* leaf_nodes,
                          IndexType* offsets,
                          IndexType* counts,
                          IndexType*& candidates,
                          IndexType num_objs,
                          const PrimitiveType* objs,
                          int allocatorID)
{
  SLIC_ASSERT(inner_nodes != nullptr);
  SLIC_ASSERT(inner_node_children != nullptr);
  SLIC_ASSERT(leaf_nodes != nullptr);
  SLIC_ERROR_IF(offsets == nullptr, "supplied null pointer for offsets!");
  SLIC_ERROR_IF(counts == nullptr, "supplied null pointer for counts!");
  SLIC_ERROR_IF(objs == nullptr, "supplied null pointer for test primitives!");

  int total_count = 0;

  // STEP 1: count number of candidates for each query point
  using reduce_pol = typename axom::execution_space<ExecSpace>::reduce_policy;
  RAJA::ReduceSum<reduce_pol, IndexType> total_count_reduce(0);

  AXOM_PERF_MARK_SECTION("PASS[1]:count_traversal",
                         for_all<ExecSpace>(
                           N,
                           AXOM_LAMBDA(IndexType i) {
                             int32 count = 0;
                             PrimitiveType primitive = objs[i];

                             BVH_LEAF_ACTION(
                               leafAction,
                               int32 AXOM_NOT_USED(current_node),
                               const int32* AXOM_NOT_USED(leaf_nodes))
                             {
                               count++;
                             };

                             lbvh::bvh_traverse(inner_nodes,
                                                inner_node_children,
                                                leaf_nodes,
                                                primitive,
                                                predicate,
                                                leafAction);

                             counts[i] = count;
                             total_count_reduce += count;
                           }););

  // STEP 2: exclusive scan to get offsets in candidate array for each query
  using exec_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  AXOM_PERF_MARK_SECTION(
    "exclusive_scan",
    RAJA::exclusive_scan<exec_policy>(counts,
                                      counts + num_objs,
                                      offsets,
                                      RAJA::operators::plus<IndexType> {}););

  // STEP 3: allocate memory for all candidates
  AXOM_PERF_MARK_SECTION(
    "allocate_candidates",
    IndexType total_candidates = static_cast<IndexType>(total_count);
    candidates = axom::allocate<IndexType>(total_candidates, allocatorID););

  // STEP 4: fill in candidates for each point
  AXOM_PERF_MARK_SECTION(
    "PASS[2]:fill_traversal",
    for_all<ExecSpace>(
      numObjs,
      AXOM_LAMBDA(IndexType i) {
        int32 offset = offsets[i];

        PrimitiveType obj = objs[i];

        BVH_LEAF_ACTION(leafAction, int32 current_node, const int32* leaf_nodes)
        {
          candidates[offset] = leaf_nodes[current_node];
          offset++;
        };

        lbvh::bvh_traverse(inner_nodes,
                           inner_node_children,
                           leaf_nodes,
                           obj,
                           predicate,
                           leafAction);
      }););
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

  PointType* packed_points = axom::allocate<PointType>(numPts, m_AllocatorID);
  for_all<ExecSpace>(
    numPts,
    AXOM_LAMBDA(IndexType i) {
      QueryAccessor::getPoint(packed_points[i], i, x, y, z);
    });

  // STEP 2: define traversal predicates
  BVH_PREDICATE(predicate, const PointType& p, const BoundingBoxType& bb)
  {
    return bb.contains(p);
  };

  find_candidates_impl<NDIMS, ExecSpace>(predicate,
                                         inner_nodes,
                                         inner_node_children,
                                         leaf_nodes,
                                         offsets,
                                         counts,
                                         candidates,
                                         numPts,
                                         packed_points,
                                         m_AllocatorID);

  axom::deallocate(packed_points);
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

  RayType* packed_rays = axom::allocate<RayType>(numRays, m_AllocatorID);
  for_all<ExecSpace>(
    numRays,
    AXOM_LAMBDA(IndexType i) {
      typename RayType::PointType origin;
      typename RayType::VectorType direction;
      QueryAccessor::getPoint(origin, i, x0, y0, z0);
      QueryAccessor::getPoint(direction, i, nx, ny, nz);

      RayType ray {origin, direction};
      packed_rays[i] = ray;
    });

  // STEP 2: define traversal predicates
  BVH_PREDICATE(predicate, const RayType& r, const BoundingBoxType& bb)
  {
    primal::Point<FloatType, NDIMS> tmp;
    return primal::detail::intersect_ray(r, bb, tmp, TOL);
  };

  find_candidates_impl<NDIMS, ExecSpace>(predicate,
                                         inner_nodes,
                                         inner_node_children,
                                         leaf_nodes,
                                         offsets,
                                         counts,
                                         candidates,
                                         numRays,
                                         packed_rays,
                                         m_AllocatorID);

  axom::deallocate(packed_rays);
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

  BoundingBoxType* packed_boxes =
    axom::allocate<BoundingBoxType>(numBoxes, m_AllocatorID);
  for_all<ExecSpace>(
    numBoxes,
    AXOM_LAMBDA(IndexType i) {
      QueryAccessor::getBoundingBox(packed_boxes[i],
                                    i,
                                    xmin,
                                    xmax,
                                    ymin,
                                    ymax,
                                    zmin,
                                    zmax);
    });

  // STEP 2: define traversal predicates
  BVH_PREDICATE(predicate, const BoundingBoxType& bb1, const BoundingBoxType& bb2)
  {
    return bb1.intersectsWith(bb2);
  };

  find_candidates_impl<NDIMS, ExecSpace>(predicate,
                                         inner_nodes,
                                         inner_node_children,
                                         leaf_nodes,
                                         offsets,
                                         counts,
                                         candidates,
                                         numBoxes,
                                         packed_boxes,
                                         m_AllocatorID);

  axom::deallocate(packed_boxes);
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
