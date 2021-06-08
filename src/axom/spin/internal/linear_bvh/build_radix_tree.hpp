// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BUILD_RADIX_TREE_H_
#define AXOM_SPIN_BUILD_RADIX_TREE_H_

#include "axom/config.hpp"  // for axom compile-time definitions

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"

#include "axom/core/utilities/AnnotationMacros.hpp"  // for annotations

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"

#include "axom/spin/internal/linear_bvh/BVHData.hpp"
#include "axom/spin/internal/linear_bvh/RadixTree.hpp"

#include "axom/spin/MortonIndex.hpp"

#include "axom/core/utilities/Utilities.hpp"  // for isNearlyEqual()
#include "axom/slic/interface/slic.hpp"       // for slic

// RAJA includes
#include "RAJA/RAJA.hpp"

#if defined(AXOM_USE_CUDA)
  // NOTE: uses the cub installation that is  bundled with RAJA
  #include "cub/device/device_radix_sort.cuh"
#endif

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{
//------------------------------------------------------------------------------
//Returns 30 bit morton code for coordinates for
// x, y, and z are expecting to be between [0,1]
template <typename FloatType, int Dims>
static inline AXOM_HOST_DEVICE axom::int32 morton32_encode(
  const primal::Vector<FloatType, Dims>& point)
{
  //for a float, take the first 10 bits. Note, 2^10 = 1024
  constexpr int NUM_BITS_PER_DIM = 32 / Dims;
  constexpr FloatType FLOAT_TO_INT = 1 << NUM_BITS_PER_DIM;
  constexpr FloatType FLOAT_CEILING = FLOAT_TO_INT - 1;

  int32 int_coords[Dims];
  for(int i = 0; i < Dims; i++)
  {
    int_coords[i] =
      fmin(fmax(point[i] * FLOAT_TO_INT, (FloatType)0), FLOAT_CEILING);
  }

  primal::Point<int32, Dims> integer_pt(int_coords);

  return convertPointToMorton<int32>(integer_pt);
}

//------------------------------------------------------------------------------
//Returns 30 bit morton code for coordinates for
//coordinates in the unit cude
static inline AXOM_HOST_DEVICE axom::int64 morton64_encode(axom::float32 x,
                                                           axom::float32 y,
                                                           axom::float32 z = 0.0)
{
  //take the first 21 bits. Note, 2^21= 2097152.0f
  x = fmin(fmax(x * 2097152.0f, 0.0f), 2097151.0f);
  y = fmin(fmax(y * 2097152.0f, 0.0f), 2097151.0f);
  z = fmin(fmax(z * 2097152.0f, 0.0f), 2097151.0f);

  primal::Point<int64, 3> integer_pt =
    primal::Point<int64, 3>::make_point((int64)x, (int64)y, (int64)z);

  return convertPointToMorton<int64>(integer_pt);
}

template <typename ExecSpace, typename FloatType, int NDIMS>
void transform_boxes(const FloatType* boxes,
                     primal::BoundingBox<FloatType, NDIMS>* aabbs,
                     int32 size,
                     FloatType scale_factor)
{
  AXOM_PERF_MARK_FUNCTION("transform_boxes");

  constexpr int STRIDE = 2 * NDIMS;

  for_all<ExecSpace>(
    size,
    AXOM_LAMBDA(int32 i) {
      primal::BoundingBox<FloatType, NDIMS> aabb;

      const int32 offset = i * STRIDE;

      primal::Point<FloatType, NDIMS> min_point(boxes + offset);
      primal::Point<FloatType, NDIMS> max_point(boxes + offset + NDIMS);

      aabb.addPoint(min_point);
      aabb.addPoint(max_point);
      aabb.scale(scale_factor);

      aabbs[i] = aabb;
    });
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType, int NDIMS>
primal::BoundingBox<FloatType, NDIMS> reduce(
  primal::BoundingBox<FloatType, NDIMS>* aabbs,
  int32 size)
{
  AXOM_PERF_MARK_FUNCTION("reduce_abbs");

  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

  primal::Point<FloatType, NDIMS> min_pt, max_pt;

  FloatType infinity = std::numeric_limits<FloatType>::max();
  FloatType neg_infinity = std::numeric_limits<FloatType>::lowest();

  for(int dim = 0; dim < NDIMS; dim++)
  {
    RAJA::ReduceMin<reduce_policy, FloatType> min_coord(infinity);
    RAJA::ReduceMax<reduce_policy, FloatType> max_coord(neg_infinity);

    for_all<ExecSpace>(
      size,
      AXOM_LAMBDA(int32 i) {
        const primal::BoundingBox<FloatType, NDIMS>& aabb = aabbs[i];
        min_coord.min(aabb.getMin()[dim]);
        max_coord.max(aabb.getMax()[dim]);
      });

    min_pt[dim] = min_coord.get();
    max_pt[dim] = max_coord.get();
  }

  return primal::BoundingBox<FloatType, NDIMS>(min_pt, max_pt);
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType, int NDIMS>
void get_mcodes(primal::BoundingBox<FloatType, NDIMS>* aabbs,
                int32 size,
                const primal::BoundingBox<FloatType, NDIMS>& bounds,
                uint32* mcodes)
{
  AXOM_PERF_MARK_FUNCTION("get_mcodes");

  primal::Vector<FloatType, NDIMS> extent, inv_extent, min_coord;

  extent = bounds.getMax();
  extent -= bounds.getMin();
  min_coord = bounds.getMin();

  for(int i = 0; i < NDIMS; ++i)
  {
    inv_extent[i] = utilities::isNearlyEqual<FloatType>(extent[i], .0f)
      ? 0.f
      : 1.f / extent[i];
  }

  for_all<ExecSpace>(
    size,
    AXOM_LAMBDA(int32 i) {
      const primal::BoundingBox<FloatType, NDIMS>& aabb = aabbs[i];

      // get the center and normalize it
      primal::Vector<FloatType, NDIMS> centroid = aabb.getCentroid();
      centroid = (centroid - min_coord).array() * inv_extent.array();
      mcodes[i] = morton32_encode(centroid);
    });
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename IntType>
void array_counting(IntType* iterator,
                    const IntType& size,
                    const IntType& start,
                    const IntType& step)
{
  AXOM_PERF_MARK_FUNCTION("array_counting");

  for_all<ExecSpace>(
    size,
    AXOM_LAMBDA(int32 i) { iterator[i] = start + i * step; });
}

//------------------------------------------------------------------------------
//
// reorder and array based on a new set of indices.
// array   [a,b,c]
// indices [1,0,2]
// result  [b,a,c]
//
template <typename ExecSpace, typename T>
void reorder(int32* indices, T*& array, int32 size, int allocatorID)
{
  AXOM_PERF_MARK_FUNCTION("reorder");

  T* temp = axom::allocate<T>(size, allocatorID);

  for_all<ExecSpace>(
    size,
    AXOM_LAMBDA(int32 i) {
      int32 in_idx = indices[i];
      temp[i] = array[in_idx];
    });

  axom::deallocate(array);
  array = temp;
}

//------------------------------------------------------------------------------
#if(RAJA_VERSION_MAJOR > 0) || \
  ((RAJA_VERSION_MAJOR == 0) && (RAJA_VERSION_MINOR >= 12))

template <typename ExecSpace>
void sort_mcodes(uint32*& mcodes, int32 size, int32* iter)
{
  AXOM_PERF_MARK_FUNCTION("sort_mcodes");

  array_counting<ExecSpace>(iter, size, 0, 1);

  AXOM_PERF_MARK_SECTION(
    "raja_stable_sort",
    using EXEC_POL = typename axom::execution_space<ExecSpace>::loop_policy;
    RAJA::stable_sort_pairs<EXEC_POL>(mcodes, mcodes + size, iter););
}

#else

// fall back to std::stable_sort
template <typename ExecSpace>
void sort_mcodes(uint32*& mcodes, int32 size, int32* iter)
{
  AXOM_PERF_MARK_FUNCTION("sort_mcodes");

  array_counting<ExecSpace>(iter, size, 0, 1);

  AXOM_PERF_MARK_SECTION(
    "cpu_sort",

    std::stable_sort(iter,
                     iter + size,
                     [=](int32 i1, int32 i2) { return mcodes[i1] < mcodes[i2]; });

  );

  const int allocID = axom::execution_space<ExecSpace>::allocatorID();
  reorder<ExecSpace>(iter, mcodes, size, allocID);
}

#endif /* RAJA version 0.12.0 and above */

//------------------------------------------------------------------------------
//
// count leading zeros
//
inline AXOM_HOST_DEVICE axom::int32 clz(axom::int32 x)
{
  axom::int32 y;
  axom::int32 n = 32;
  y = x >> 16;
  if(y != 0)
  {
    n = n - 16;
    x = y;
  }
  y = x >> 8;
  if(y != 0)
  {
    n = n - 8;
    x = y;
  }
  y = x >> 4;
  if(y != 0)
  {
    n = n - 4;
    x = y;
  }
  y = x >> 2;
  if(y != 0)
  {
    n = n - 2;
    x = y;
  }
  y = x >> 1;
  if(y != 0) return axom::int32(n - 2);
  return axom::int32(n - x);
}

//------------------------------------------------------------------------------
template <typename IntType, typename MCType>
AXOM_HOST_DEVICE IntType delta(const IntType& a,
                               const IntType& b,
                               const IntType& inner_size,
                               const MCType* mcodes)
{
  bool tie = false;
  bool out_of_range = (b < 0 || b > inner_size);
  //still make the call but with a valid adderss
  const int32 bb = (out_of_range) ? 0 : b;
  const uint32 acode = mcodes[a];
  const uint32 bcode = mcodes[bb];
  //use xor to find where they differ
  uint32 exor = acode ^ bcode;
  tie = (exor == 0);
  //break the tie, a and b must always differ
  exor = tie ? uint32(a) ^ uint32(bb) : exor;
  int32 count = clz(exor);
  if(tie) count += 32;
  count = (out_of_range) ? -1 : count;
  return count;
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType, int NDIMS>
void build_tree(RadixTree<FloatType, NDIMS>& data)
{
  AXOM_PERF_MARK_FUNCTION("build_tree");

  // http://research.nvidia.com/sites/default/files/publications/karras2012hpg_paper.pdf

  // Pointers and vars are redeclared because I have a faint memory
  // of a huge amount of pain and suffering due so cuda
  // lambda captures of pointers inside a struct. Bad memories
  // of random segfaults ........ be warned
  const int32 inner_size = data.m_inner_size;
  int32* lchildren_ptr = data.m_left_children;
  int32* rchildren_ptr = data.m_right_children;
  int32* parent_ptr = data.m_parents;
  const uint32* mcodes_ptr = data.m_mcodes;

  for_all<ExecSpace>(
    inner_size,
    AXOM_LAMBDA(int32 i) {
      //determine range direction
      int32 d = 0 > (delta(i, i + 1, inner_size, mcodes_ptr) -
                     delta(i, i - 1, inner_size, mcodes_ptr))
        ? -1
        : 1;

      //find upper bound for the length of the range
      int32 min_delta = delta(i, i - d, inner_size, mcodes_ptr);
      int32 lmax = 2;
      while(delta(i, i + lmax * d, inner_size, mcodes_ptr) > min_delta)
        lmax *= 2;

      //binary search to find the lower bound
      int32 l = 0;
      for(int32 t = lmax / 2; t >= 1; t /= 2)
      {
        if(delta(i, i + (l + t) * d, inner_size, mcodes_ptr) > min_delta)
        {
          l += t;
        }
      }

      int32 j = i + l * d;
      int32 delta_node = delta(i, j, inner_size, mcodes_ptr);
      int32 s = 0;
      FloatType div_factor = 2.f;
      //find the split postition using a binary search
      for(int32 t = (int32)ceil(float32(l) / div_factor);;
          div_factor *= 2, t = (int32)ceil(float32(l) / div_factor))
      {
        if(delta(i, i + (s + t) * d, inner_size, mcodes_ptr) > delta_node)
        {
          s += t;
        }
        if(t == 1) break;
      }

      int32 split = i + s * d + utilities::min(d, 0);
      // assign parent/child pointers
      if(utilities::min(i, j) == split)
      {
        //leaf
        parent_ptr[split + inner_size] = i;
        lchildren_ptr[i] = split + inner_size;
      }
      else
      {
        //inner node
        parent_ptr[split] = i;
        lchildren_ptr[i] = split;
      }

      if(utilities::max(i, j) == split + 1)
      {
        //leaf
        parent_ptr[split + inner_size + 1] = i;
        rchildren_ptr[i] = split + inner_size + 1;
      }
      else
      {
        parent_ptr[split + 1] = i;
        rchildren_ptr[i] = split + 1;
      }

      if(i == 0)
      {
        // flag the root
        parent_ptr[0] = -1;
      }
    });
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename T>
static void array_memset(T* array, const int32 size, const T val)
{
  AXOM_PERF_MARK_FUNCTION("array_memset");

  for_all<ExecSpace>(
    size,
    AXOM_LAMBDA(int32 i) { array[i] = val; });
}

/*!
 * \def SPIN_BVH_THREAD_FENCE_SYSTEM
 *
 * \brief Macro for __threadfence_system() on NVIDIA GPUs expands
 *  to nothing on all other platforms.
 *
 * \note This is an internal macro use in the propagate_aabbs() method below.
 */
#if defined(__CUDA_ARCH__)
  // on NVIDIA GPUs call __threadfence_system()
  #define SPIN_BVH_THREAD_FENCE_SYSTEM __threadfence_system()
#else
  // on CPUs do nothing
  #define SPIN_BVH_THREAD_FENCE_SYSTEM
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType, int NDIMS>
void propagate_aabbs(RadixTree<FloatType, NDIMS>& data, int allocatorID)
{
  AXOM_PERF_MARK_FUNCTION("propagate_abbs");

  const int inner_size = data.m_inner_size;
  const int leaf_size = data.m_inner_size + 1;
  SLIC_ASSERT(leaf_size == data.m_size);

  // Pointers and vars are redeclared because I have a faint memory
  // of a huge amount of pain and suffering due so cuda
  // labda captures of pointers indide a struct. Bad memories
  // of random segfaults ........ be warned
  const int32* lchildren_ptr = data.m_left_children;
  const int32* rchildren_ptr = data.m_right_children;
  const int32* parent_ptr = data.m_parents;
  const primal::BoundingBox<FloatType, NDIMS>* leaf_aabb_ptr = data.m_leaf_aabbs;

  primal::BoundingBox<FloatType, NDIMS>* inner_aabb_ptr = data.m_inner_aabbs;

  int32* counters_ptr = axom::allocate<int32>(inner_size, allocatorID);

  array_memset<ExecSpace>(counters_ptr, inner_size, 0);

  using atomic_policy = typename axom::execution_space<ExecSpace>::atomic_policy;

  for_all<ExecSpace>(
    leaf_size,
    AXOM_LAMBDA(int32 i) {
      int32 current_node = parent_ptr[inner_size + i];

      while(current_node != -1)
      {
        int32 old =
          RAJA::atomicAdd<atomic_policy>(&(counters_ptr[current_node]), 1);

        if(old == 0)
        {
          // first thread to get here kills itself
          return;
        }

        int32 lchild = lchildren_ptr[current_node];
        int32 rchild = rchildren_ptr[current_node];

        // gather the aabbs
        primal::BoundingBox<FloatType, NDIMS> aabb;
        if(lchild >= inner_size)
        {
          aabb.addBox(leaf_aabb_ptr[lchild - inner_size]);
        }
        else
        {
          aabb.addBox(inner_aabb_ptr[lchild]);
        }

        if(rchild >= inner_size)
        {
          aabb.addBox(leaf_aabb_ptr[rchild - inner_size]);
        }
        else
        {
          aabb.addBox(inner_aabb_ptr[rchild]);
        }

        inner_aabb_ptr[current_node] = aabb;

        // NOTE: this ensures that the write to global memory above
        // is observed by all other threads. It provides a workaround
        // where a thread would write to the L1 cache instead, which
        // would not be observed by other threads.
        //
        // See Github issue: https://github.com/LLNL/axom/issues/307
        //
        // This is a workaround for NVIDIA GPUs.
        SPIN_BVH_THREAD_FENCE_SYSTEM;

        current_node = parent_ptr[current_node];
      }
    });

  axom::deallocate(counters_ptr);
}

#undef SPIN_BVH_THREAD_FENCE_SYSTEM

//------------------------------------------------------------------------------
template <typename ExecSpace, typename FloatType, int NDIMS>
void build_radix_tree(const FloatType* boxes,
                      int size,
                      primal::BoundingBox<FloatType, NDIMS>& bounds,
                      RadixTree<FloatType, NDIMS>& radix_tree,
                      FloatType scale_factor,
                      int allocatorID)
{
  AXOM_PERF_MARK_FUNCTION("build_radix_tree");

  // sanity checks
  SLIC_ASSERT(boxes != nullptr);
  SLIC_ASSERT(size > 0);

  radix_tree.allocate(size, allocatorID);

  // copy so we don't reorder the input
  transform_boxes<ExecSpace>(boxes, radix_tree.m_leaf_aabbs, size, scale_factor);

  // evaluate global bounds
  bounds = reduce<ExecSpace>(radix_tree.m_leaf_aabbs, size);

  // sort aabbs based on morton code
  // original positions of the sorted morton codes.
  // allows us to gather / sort other arrays.
  get_mcodes<ExecSpace>(radix_tree.m_leaf_aabbs, size, bounds, radix_tree.m_mcodes);
  sort_mcodes<ExecSpace>(radix_tree.m_mcodes, size, radix_tree.m_leafs);

  reorder<ExecSpace>(radix_tree.m_leafs, radix_tree.m_leaf_aabbs, size, allocatorID);

  build_tree<ExecSpace>(radix_tree);

  propagate_aabbs<ExecSpace>(radix_tree, allocatorID);
}

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */
#endif
