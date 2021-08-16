// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_POLICY_LINEARBVH_HPP_
#define AXOM_SPIN_POLICY_LINEARBVH_HPP_

// axom core includes
#include "axom/core/Types.hpp"              // for fixed bitwidth types
#include "axom/core/execution/for_all.hpp"  // for generic for_all()
#include "axom/core/memory_management.hpp"  // for alloc()/free()

#include "axom/core/utilities/AnnotationMacros.hpp"  // for annotations

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Vector.hpp"

// linear bvh includes
#include "axom/spin/internal/linear_bvh/RadixTree.hpp"
#include "axom/spin/internal/linear_bvh/build_radix_tree.hpp"
#include "axom/spin/internal/linear_bvh/bvh_traverse.hpp"
#include "axom/spin/internal/linear_bvh/bvh_vtkio.hpp"

// C/C++ includes
#include <fstream>  // for std::ofstream
#include <sstream>  // for std::ostringstream
#include <string>   // for std::string
#include <vector>   // for std::vector

namespace axom
{
namespace spin
{
namespace policy
{
namespace lbvh = internal::linear_bvh;

/*!
 * \brief LinearBVH provides a policy for a BVH implementation which supports
 *  parallel linear construction on both CPU and GPU.
 *
 * \note The internal node data is organized into two arrays, one containing
 *  bounding boxes of the internal nodes' two child nodes, and the other
 *  containing indices of the two child nodes (scaled by two if inner node,
 *  ones-complement if leaf node).
 *
 */
template <typename FloatType, int NDIMS, typename ExecSpace>
class LinearBVH
{
public:
  using BoundingBoxType = primal::BoundingBox<FloatType, NDIMS>;

  LinearBVH() = default;
  ~LinearBVH() { deallocate(); }

  /*!
   * \brief Builds a linear BVH with the given bounding boxes as leaf nodes.
   *
   * \param [in] boxes the bounding boxes for each leaf node
   * \param [in] numBoxes the number of bounding boxes
   * \param [in] scaleFactor scale factor applied to each bounding box before insertion into the BVH
   */
  void buildImpl(const BoundingBoxType* boxes,
                 IndexType numBoxes,
                 FloatType scaleFactor,
                 int allocatorID);

  /*!
   * \brief Performs a traversal to find the candidates for each query primitive.
   *
   * \param [in] predicate traversal predicate functor for bin check.
   * \param [out] offsets array of offsets into the candidate array for each query primitive
   * \param [out] counts array of candidate counts for each query primitive
   * \param [out] candidates array of the potential candidates for intersection with the BVH
   * \param [in] numObjs the number of user-supplied query primitives
   * \param [in] objs array of primitives to query against the BVH
   *
   * \return total_count the total count of candidates for all query primitives.
   */
  template <typename PrimitiveType, typename Predicate, typename PrimitiveIndexable>
  void findCandidatesImpl(Predicate&& predicate,
                          IndexType* offsets,
                          IndexType* counts,
                          IndexType*& candidates,
                          IndexType numObjs,
                          PrimitiveIndexable objs,
                          int allocatorID) const;

  void writeVtkFileImpl(const std::string& fileName) const;

  BoundingBoxType getBoundsImpl() const { return m_bounds; }

private:
  void allocate(int32 size, int allocID)
  {
    AXOM_PERF_MARK_FUNCTION("LinearBVH::allocate");
    m_inner_nodes = axom::allocate<BoundingBoxType>((size - 1) * 2, allocID);
    m_inner_node_children = axom::allocate<int32>((size - 1) * 2, allocID);
    m_leaf_nodes = axom::allocate<int32>(size, allocID);
  }

  void deallocate()
  {
    AXOM_PERF_MARK_FUNCTION("LinearBVH::deallocate");
    axom::deallocate(m_inner_nodes);
    axom::deallocate(m_inner_node_children);
    axom::deallocate(m_leaf_nodes);
  }

  BoundingBoxType* m_inner_nodes {nullptr};  // BVH bins including leafs
  int32* m_inner_node_children {nullptr};
  int32* m_leaf_nodes {nullptr};  // leaf data
  primal::BoundingBox<FloatType, NDIMS> m_bounds;
};

template <typename FloatType, int NDIMS, typename ExecSpace>
void LinearBVH<FloatType, NDIMS, ExecSpace>::buildImpl(const BoundingBoxType* boxes,
                                                       IndexType numBoxes,
                                                       FloatType scaleFactor,
                                                       int allocatorID)
{
  AXOM_PERF_MARK_FUNCTION("LinearBVH::buildImpl");

  // STEP 1: Build a RadixTree consisting of the bounding boxes, sorted
  // by their corresponding morton code.
  lbvh::RadixTree<FloatType, NDIMS> radix_tree;
  primal::BoundingBox<FloatType, NDIMS> global_bounds;
  lbvh::build_radix_tree<ExecSpace>(boxes,
                                    numBoxes,
                                    global_bounds,
                                    radix_tree,
                                    scaleFactor,
                                    allocatorID);

  // STEP 2: emit the BVH data-structure from the radix tree
  m_bounds = global_bounds;
  allocate(numBoxes, allocatorID);

  // STEP 3: emit the BVH
  const int32 size = radix_tree.m_size;
  const int32 inner_size = radix_tree.m_inner_size;
  SLIC_ASSERT(inner_size == size - 1);

  const int32* lchildren_ptr = radix_tree.m_left_children;
  const int32* rchildren_ptr = radix_tree.m_right_children;

  const BoundingBoxType* leaf_aabb_ptr = radix_tree.m_leaf_aabbs;
  const BoundingBoxType* inner_aabb_ptr = radix_tree.m_inner_aabbs;

  primal::BoundingBox<FloatType, NDIMS>* bvh_inner_nodes = m_inner_nodes;
  int32* bvh_inner_node_children = m_inner_node_children;

  AXOM_PERF_MARK_SECTION("emit_bvh_parents",
                         for_all<ExecSpace>(
                           inner_size,
                           AXOM_LAMBDA(int32 node) {
                             BoundingBoxType l_aabb, r_aabb;

                             int32 lchild = lchildren_ptr[node];
                             if(lchild >= inner_size)
                             {
                               l_aabb = leaf_aabb_ptr[lchild - inner_size];
                               lchild = -(lchild - inner_size + 1);
                             }
                             else
                             {
                               l_aabb = inner_aabb_ptr[lchild];
                               // do the offset now
                               lchild *= 2;
                             }

                             int32 rchild = rchildren_ptr[node];
                             if(rchild >= inner_size)
                             {
                               r_aabb = leaf_aabb_ptr[rchild - inner_size];
                               rchild = -(rchild - inner_size + 1);
                             }
                             else
                             {
                               r_aabb = inner_aabb_ptr[rchild];
                               // do the offset now
                               rchild *= 2;
                             }

                             const int32 out_offset = node * 2;
                             bvh_inner_nodes[out_offset + 0] = l_aabb;
                             bvh_inner_nodes[out_offset + 1] = r_aabb;

                             bvh_inner_node_children[out_offset + 0] = lchild;
                             bvh_inner_node_children[out_offset + 1] = rchild;
                           }););

  int32* radix_tree_leafs = radix_tree.m_leafs;
  int32* bvh_leafs = m_leaf_nodes;

  AXOM_PERF_MARK_SECTION(
    "emit_bvh_leafs",
    for_all<ExecSpace>(
      size,
      AXOM_LAMBDA(int32 i) { bvh_leafs[i] = radix_tree_leafs[i]; }););

  radix_tree.deallocate();
}

template <typename FloatType, int NDIMS, typename ExecSpace>
template <typename PrimitiveType, typename Predicate, typename PrimitiveIndexable>
void LinearBVH<FloatType, NDIMS, ExecSpace>::findCandidatesImpl(
  Predicate&& predicate,
  IndexType* offsets,
  IndexType* counts,
  IndexType*& candidates,
  IndexType numObjs,
  PrimitiveIndexable objs,
  int allocatorID) const

{
  AXOM_PERF_MARK_FUNCTION("LinearBVH::findCandidatesImpl");

  SLIC_ERROR_IF(offsets == nullptr, "supplied null pointer for offsets!");
  SLIC_ERROR_IF(counts == nullptr, "supplied null value for counts!");

  const BoundingBoxType* inner_nodes = m_inner_nodes;
  const int32* inner_node_children = m_inner_node_children;
  const int32* leaf_nodes = m_leaf_nodes;
  SLIC_ASSERT(inner_nodes != nullptr);
  SLIC_ASSERT(inner_node_children != nullptr);
  SLIC_ASSERT(leaf_nodes != nullptr);

#if defined(AXOM_USE_RAJA)
  // STEP 1: count number of candidates for each query point
  using reduce_pol = typename axom::execution_space<ExecSpace>::reduce_policy;
  RAJA::ReduceSum<reduce_pol, IndexType> total_count_reduce(0);

  AXOM_PERF_MARK_SECTION(
    "PASS[1]:count_traversal",
    for_all<ExecSpace>(
      numObjs,
      AXOM_LAMBDA(IndexType i) {
        int32 count = 0;
        PrimitiveType primitive {objs[i]};

        auto leafAction = [&count](int32 AXOM_NOT_USED(current_node),
                                   const int32* AXOM_NOT_USED(leaf_nodes)) {
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
                                      counts + numObjs,
                                      offsets,
                                      RAJA::operators::plus<IndexType> {}););

  int total_count = total_count_reduce.get();

  // STEP 3: allocate memory for all candidates
  AXOM_PERF_MARK_SECTION(
    "allocate_candidates",
    IndexType total_candidates = static_cast<IndexType>(total_count);
    candidates = axom::allocate<IndexType>(total_candidates, allocatorID););

  // STEP 4: fill in candidates for each point
  AXOM_PERF_MARK_SECTION("PASS[2]:fill_traversal",
                         for_all<ExecSpace>(
                           numObjs,
                           AXOM_LAMBDA(IndexType i) {
                             int32 offset = offsets[i];

                             PrimitiveType obj {objs[i]};
                             auto leafAction = [&offset, candidates](
                                                 int32 current_node,
                                                 const int32* leafs) {
                               candidates[offset] = leafs[current_node];
                               offset++;
                             };

                             lbvh::bvh_traverse(inner_nodes,
                                                inner_node_children,
                                                leaf_nodes,
                                                obj,
                                                predicate,
                                                leafAction);
                           }););
#else  // CPU-only and no RAJA: do single traversal

  std::vector<int> search_candidates;
  int current_offset = 0;

  // STEP 1: do single-pass traversal with std::vector for candidates
  AXOM_PERF_MARK_SECTION(
    "PASS[1]:fill_traversal", for_all<ExecSpace>(numObjs, [&](IndexType i) {
      int matching_leaves = 0;
      PrimitiveType obj {objs[i]};
      offsets[i] = current_offset;

      auto leafAction = [&](int32 current_node, const int32* leafs) {
        search_candidates.emplace_back(leafs[current_node]);
        matching_leaves++;
        current_offset++;
      };

      lbvh::bvh_traverse(inner_nodes,
                         inner_node_children,
                         leaf_nodes,
                         obj,
                         predicate,
                         leafAction);
      counts[i] = matching_leaves;
    }););

  SLIC_ASSERT(current_offset == static_cast<IndexType>(search_candidates.size()));

  // STEP 2: allocate memory for all candidates
  AXOM_PERF_MARK_SECTION(
    "allocate_candidates", IndexType total_candidates = search_candidates.size();
    candidates = axom::allocate<IndexType>(total_candidates, allocatorID););

  // STEP 3: copy candiates to destination array
  std::copy(search_candidates.begin(), search_candidates.end(), candidates);
#endif
}

template <typename FloatType, int NDIMS, typename ExecSpace>
void LinearBVH<FloatType, NDIMS, ExecSpace>::writeVtkFileImpl(
  const std::string& fileName) const
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
  lbvh::write_root(m_bounds, numPoints, numBins, nodes, cells, levels);

  // STEP 2: traverse the BVH and dump each bin
  constexpr int32 ROOT = 0;
  lbvh::write_recursive<FloatType, NDIMS>(m_inner_nodes,
                                          m_inner_node_children,
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

}  // namespace policy
}  // namespace spin
}  // namespace axom
#endif /* AXOM_SPIN_POLICY_LINEARBVH_HPP_ */
