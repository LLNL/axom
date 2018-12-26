#ifndef AXOM_PRIMAL_BVH_LINEAR_BUILDER_H_
#define AXOM_PRIMAL_BVH_LINEAR_BUILDER_H_

#include "axom/core/Types.hpp"    // for fixed bitwidth types
#include "axom/primal/spatial_acceleration/linear_bvh/aabb.hpp"

#include "umpire/Umpire.hpp"
using namespace axom::common;

namespace axom
{
namespace primal
{
namespace bvh
{
static std::string alloc_name = "HOST";

template<typename T>
T * umpire_alloc(int64 size)
{
 auto& rm = umpire::ResourceManager::getInstance();
 umpire::Allocator allocator = rm.getAllocator(alloc_name);
 return static_cast<T*>(allocator.allocate(size*sizeof(T)));
}

template<typename T> void umpire_dealloc(T* ptr)
{
  auto& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator(alloc_name);
  allocator.deallocate(ptr);
}

template<typename T> void umpire_copy(T *src, T *dest)
{
  auto& rm = umpire::ResourceManager::getInstance();
  rm.copy(src, dest);
}

struct BVH
{
  Vec<float32, 4>  *m_inner_nodes;
  int32            *m_leaf_nodes;
  AABB              m_bounds;

  BVH()
  {
    m_inner_nodes = nullptr;
    m_leaf_nodes = nullptr;
  }

  ~BVH()
  {
    umpire_dealloc(m_inner_nodes);
    umpire_dealloc(m_leaf_nodes);
  }
};

class LinearBVHBuilder
{

public:
  BVH construct(const double *boxes, int size);

};

} /* namespace axom */
} /* namespace primal */
} /* namespace bvh */
#endif
