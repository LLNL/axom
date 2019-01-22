/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "axom/primal/spatial_acceleration/linear_bvh/linear_bvh_builder.hpp"

#include "axom/primal/spatial_acceleration/linear_bvh/math.hpp"
#include "axom/primal/spatial_acceleration/linear_bvh/morton_codes.hpp"
#include "axom/primal/spatial_acceleration/linear_bvh/policies.hpp"

// RAJA includes
#include "RAJA/RAJA.hpp"

using namespace axom::common;

namespace axom
{
namespace primal
{
namespace bvh
{

void transform_boxes(const double *boxes, AABB *aabbs, const int32 size)
{

  RAJA::forall<raja_for_policy>(RAJA::RangeSegment(0, size), AXOM_LAMBDA (int32 i)
  {
    AABB aabb;
    Vec<float32,3> min_point, max_point;
    const int32 offset = i * 6;
    min_point[0] = boxes[offset + 0];
    min_point[1] = boxes[offset + 1];
    min_point[2] = boxes[offset + 2];
    max_point[0] = boxes[offset + 3];
    max_point[1] = boxes[offset + 4];
    max_point[2] = boxes[offset + 5];
    aabb.include(min_point);
    aabb.include(max_point);
    aabbs[i] = aabb;
  });

}

AABB reduce(AABB *aabbs, const int size)
{

  RAJA::ReduceMin<raja_reduce_policy, float32> xmin(infinity32());
  RAJA::ReduceMin<raja_reduce_policy, float32> ymin(infinity32());
  RAJA::ReduceMin<raja_reduce_policy, float32> zmin(infinity32());

  RAJA::ReduceMax<raja_reduce_policy, float32> xmax(neg_infinity32());
  RAJA::ReduceMax<raja_reduce_policy, float32> ymax(neg_infinity32());
  RAJA::ReduceMax<raja_reduce_policy, float32> zmax(neg_infinity32());

  RAJA::forall<raja_for_policy>(RAJA::RangeSegment(0, size), AXOM_LAMBDA (int32 i)
  {

    const AABB &aabb = aabbs[i];
    //std::cout<<i<<" "<<aabb<<"\n";
    xmin.min(aabb.m_x.min());
    ymin.min(aabb.m_y.min());
    zmin.min(aabb.m_z.min());

    xmax.max(aabb.m_x.max());
    ymax.max(aabb.m_y.max());
    zmax.max(aabb.m_z.max());

  });

  AABB res;
  Vec3f mins = make_vec3f(xmin.get(), ymin.get(), zmin.get());
  Vec3f maxs = make_vec3f(xmax.get(), ymax.get(), zmax.get());

  res.include(mins);
  res.include(maxs);
  return res;
}

uint32 *get_mcodes(AABB *aabbs, const int32  size, const AABB &bounds)
{
  Vec3f extent, inv_extent, min_coord;
  extent[0] = bounds.m_x.max() - bounds.m_x.min();
  extent[1] = bounds.m_y.max() - bounds.m_y.min();
  extent[2] = bounds.m_z.max() - bounds.m_z.min();

  min_coord[0] = bounds.m_x.min();
  min_coord[1] = bounds.m_y.min();
  min_coord[2] = bounds.m_z.min();

  for(int i = 0; i < 3; ++i)
  {
    inv_extent[i] = (extent[i] == .0f) ? 0.f : 1.f / extent[i];
  }

  uint32 *mcodes = axom::alloc<uint32>(size);

  RAJA::forall<raja_for_policy>(RAJA::RangeSegment(0, size), AXOM_LAMBDA (int32 i)
  {
    const AABB &aabb = aabbs[i];
    // get the center and normalize it
    float32 centroid_x = (aabb.m_x.center() - min_coord[0]) * inv_extent[0];
    float32 centroid_y = (aabb.m_y.center() - min_coord[1]) * inv_extent[1];
    float32 centroid_z = (aabb.m_z.center() - min_coord[2]) * inv_extent[2];
    mcodes[i] = morton32_encode(centroid_x, centroid_y, centroid_z);
  });

  return mcodes;
}

int32*
array_counting(const int32 &size,
               const int32 &start,
               const int32 &step)
{

  int32 *iterator = axom::alloc<int32>(size);

  RAJA::forall<raja_for_policy>(RAJA::RangeSegment(0, size), AXOM_LAMBDA (int32 i)
  {
    iterator[i] = start + i * step;
  });

  return iterator;

}
//
// reorder and array based on a new set of indices.
// array   [a,b,c]
// indices [1,0,2]
// result  [b,a,c]
//
template<typename T>
void reorder(int32 *indices, T *&array, const int32 size)
{
  T* temp = axom::alloc<T>(size);

  RAJA::forall<raja_for_policy>(RAJA::RangeSegment(0, size), AXOM_LAMBDA (int32 i)
  {
    int32 in_idx = indices[i];
    temp[i] = array[in_idx];
  });


  axom::free(array);
  array = temp;

}

int32* sort_mcodes(uint32 *&mcodes, const int32 size)
{
  int32* iter = array_counting(size, 0, 1);
  // TODO: create custom sort for GPU / CPU
  // WARNING: this will segfault with CUDA and not unified memory
  std::sort(iter,
            iter + size,
            [=](int32 i1, int32 i2)
            {
              return mcodes[i1] < mcodes[i2];
            });


  reorder(iter, mcodes, size);

  return iter;
}

struct BVHData
{
  int32  m_inner_size;
  int32  *m_left_children;
  int32  *m_right_children;
  int32  *m_parents;
  int32  *m_leafs;
  uint32 *m_mcodes;
  AABB   *m_inner_aabbs;
  AABB   *m_leaf_aabbs;
  void
  allocate(const int32 size)
  {
    m_inner_size = size;
    m_left_children = axom::alloc<int32>(size);
    m_right_children = axom::alloc<int32>(size);
    m_parents = axom::alloc<int32>(size);
    m_inner_aabbs = axom::alloc<AABB>(size);
  }

  void
  deallocate()
  {
    axom::free(m_left_children);
    axom::free(m_right_children);
    axom::free(m_parents);
    axom::free(m_inner_aabbs);
  }
};


AXOM_HOST_DEVICE int32 delta(const int32 &a,
                             const int32 &b,
                             const int32 &inner_size,
                             const uint32 *mcodes)
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
  if (tie)
    count += 32;
  count = (out_of_range) ? - 1 : count;
  return count;
}

void build_tree(BVHData &data)
{
  // http://research.nvidia.com/sites/default/files/publications/karras2012hpg_paper.pdf

  // Pointers and vars are redeclared because I have a faint memory
  // of a huge amount of pain and suffering due so cuda
  // labda captures of pointers indide a struct. Bad memories
  // of random segfaults ........ be warned
  const int32 inner_size = data.m_inner_size;
  int32 *lchildren_ptr = data.m_left_children;
  int32 *rchildren_ptr = data.m_right_children;
  int32 *parent_ptr = data.m_parents;
  const uint32 *mcodes_ptr = data.m_mcodes;

  RAJA::forall<raja_for_policy>(RAJA::RangeSegment(0, inner_size), AXOM_LAMBDA (int32 i)
  {
    //determine range direction
    int32 d = 0 > (delta(i, i + 1, inner_size, mcodes_ptr) - delta(i, i - 1, inner_size, mcodes_ptr)) ? -1 : 1;

    //find upper bound for the length of the range
    int32 min_delta = delta(i, i - d, inner_size, mcodes_ptr);
    int32 lmax = 2;
    while (delta(i, i + lmax * d, inner_size, mcodes_ptr) > min_delta)
          lmax *= 2;

    //binary search to find the lower bound
    int32 l = 0;
    for (int32 t = lmax / 2; t >= 1; t /= 2)
    {
      if (delta(i, i + (l + t) * d, inner_size, mcodes_ptr) > min_delta)
      {
        l += t;
      }
    }

    int32 j = i + l * d;
    int32 delta_node = delta(i, j, inner_size, mcodes_ptr);
    int32 s = 0;
    float32 div_factor = 2.f;
    //find the split postition using a binary search
    for (int32 t = (int32)ceil(float32(l) / div_factor);;
         div_factor *= 2, t = (int32)ceil(float32(l) / div_factor))
    {
      if (delta(i, i + (s + t) * d, inner_size, mcodes_ptr) > delta_node)
      {
        s += t;
      }
      if (t == 1) break;
    }

    int32 split = i + s * d + min(d, 0);
    // assign parent/child pointers
    if (min(i, j) == split)
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


    if (max(i, j) == split + 1)
    {
      //leaf
      parent_ptr[split + inner_size + 1] =  i;
      rchildren_ptr[i] =  split + inner_size + 1;
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

template<typename T>
static void array_memset(T* array, const int32 size, const T val)
{
  RAJA::forall<raja_for_policy>(RAJA::RangeSegment(0, size), AXOM_LAMBDA (int32 i)
  {
    array[i] = val;
  });
}

void propagate_aabbs(BVHData &data)
{
  const int inner_size = data.m_inner_size;
  const int leaf_size = data.m_inner_size + 1;
  // Pointers and vars are redeclared because I have a faint memory
  // of a huge amount of pain and suffering due so cuda
  // labda captures of pointers indide a struct. Bad memories
  // of random segfaults ........ be warned
  const int32 *lchildren_ptr = data.m_left_children;
  const int32 *rchildren_ptr = data.m_right_children;
  const int32 *parent_ptr = data.m_parents;
  const AABB  *leaf_aabb_ptr = data.m_leaf_aabbs;

  AABB  *inner_aabb_ptr = data.m_inner_aabbs;

  int32* counters_ptr = axom::alloc<int32>(inner_size);

  array_memset(counters_ptr, inner_size, 0);

  RAJA::forall<raja_for_policy>(RAJA::RangeSegment(0, leaf_size), AXOM_LAMBDA (int32 i)
  {
    int32 current_node = parent_ptr[inner_size + i];

    while(current_node != -1)
    {
      int32 old = RAJA::atomic::atomicAdd<raja_atomic_policy>(&(counters_ptr[current_node]), 1);

      if(old == 0)
      {
        // first thread to get here kills itself
        return;
      }

      int32 lchild = lchildren_ptr[current_node];
      int32 rchild = rchildren_ptr[current_node];
      // gather the aabbs
      AABB aabb;
      if(lchild >= inner_size)
      {
        aabb.include(leaf_aabb_ptr[lchild - inner_size]);
      }
      else
      {
        aabb.include(inner_aabb_ptr[lchild]);
      }

      if(rchild >= inner_size)
      {
        aabb.include(leaf_aabb_ptr[rchild - inner_size]);
      }
      else
      {
        aabb.include(inner_aabb_ptr[rchild]);
      }

      inner_aabb_ptr[current_node] = aabb;

      current_node = parent_ptr[current_node];
    }

    //printf("There can be only one\n");

  });

  axom::free(counters_ptr);
  //AABB *inner = data.m_inner_aabbs.get_host_ptr();
  //std::cout<<"Root bounds "<<inner[0]<<"\n";
}

Vec<float32,4>* emit(BVHData &data)
{
  const int inner_size = data.m_inner_size;

  const int32 *lchildren_ptr = data.m_left_children;
  const int32 *rchildren_ptr = data.m_right_children;

  const AABB  *leaf_aabb_ptr  = data.m_leaf_aabbs;
  const AABB  *inner_aabb_ptr = data.m_inner_aabbs;

  Vec<float32,4> *flat_ptr = axom::alloc<Vec<float32,4>>(inner_size * 4);;

  RAJA::forall<raja_for_policy>(RAJA::RangeSegment(0, inner_size), AXOM_LAMBDA (int32 node)
  {
    Vec<float32,4> vec1;
    Vec<float32,4> vec2;
    Vec<float32,4> vec3;
    Vec<float32,4> vec4;

    AABB l_aabb, r_aabb;

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
      lchild *= 4;
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
      rchild *= 4;
    }
    vec1[0] = l_aabb.m_x.min();
    vec1[1] = l_aabb.m_y.min();
    vec1[2] = l_aabb.m_z.min();

    vec1[3] = l_aabb.m_x.max();
    vec2[0] = l_aabb.m_y.max();
    vec2[1] = l_aabb.m_z.max();

    vec2[2] = r_aabb.m_x.min();
    vec2[3] = r_aabb.m_y.min();
    vec3[0] = r_aabb.m_z.min();

    vec3[1] = r_aabb.m_x.max();
    vec3[2] = r_aabb.m_y.max();
    vec3[3] = r_aabb.m_z.max();

    const int32 out_offset = node * 4;
    flat_ptr[out_offset + 0] = vec1;
    flat_ptr[out_offset + 1] = vec2;
    flat_ptr[out_offset + 2] = vec3;

    constexpr int32 isize = sizeof(int32);
    // memcopy so we do not truncate the ints
    memcpy(&vec4[0], &lchild, isize);
    memcpy(&vec4[1], &rchild, isize);
    flat_ptr[out_offset + 3] = vec4;
  });

  return flat_ptr;
}

BVH
LinearBVHBuilder::construct(const double *boxes, int size)
{
  // copy so we don't reorder the input
  AABB *aabbs = axom::alloc<AABB>(size);
  transform_boxes(boxes, aabbs, size);

  BVH bvh;
  AABB bounds = reduce(aabbs, size);
  uint32 *mcodes = get_mcodes(aabbs, size, bounds);

  // original positions of the sorted morton codes.
  // allows us to gather / sort other arrays.
  int32 *ids = sort_mcodes(mcodes, size);

  reorder(ids, aabbs, size);

  BVHData bvh_data;
  bvh_data.allocate(size - 1);
  // the arrays that already exist
  bvh_data.m_leafs  = ids;
  bvh_data.m_mcodes = mcodes;
  bvh_data.m_leaf_aabbs = aabbs;

  // assign parent and child pointers
  build_tree(bvh_data);

  propagate_aabbs(bvh_data);


  bvh.m_inner_nodes = emit(bvh_data);
  bvh.m_leaf_nodes = bvh_data.m_leafs;
  bvh.m_bounds = bounds;
  for(int i = 0; i < 3; ++i)
  {
    std::cout<<bvh.m_inner_nodes[i]<<"\n";
  }
  axom::free(mcodes);
  axom::free(aabbs);

  return bvh;
}
} /* namespace axom */
} /* namespace primal */
} /* namespace bvh */
