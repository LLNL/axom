// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mir/blueprint_utilities.hpp"

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
#if 0
conduit::index_t ConduitAllocateThroughAxom::conduitAllocatorID = -1;
int ConduitAllocateThroughAxom::axomAllocatorID = 0;

ConduitAllocateThroughAxom::ConduitAllocateThroughAxom(int _allocatorID)
{
  // Stash the Axom allocatorID so we don't just use the default later.
  axomAllocatorID = _allocatorID;

  // Register internal functions here as Conduit allocator functions. This returns
  // the conduitAllocatorID, which is the value that we must pass to Node::set_allocator.
  if(conduitAllocatorID == -1)
  {
    conduitAllocatorID = conduit::utils::register_allocator(internal_allocate, internal_free);
  }
}

ConduitAllocateThroughAxom::~ConduitAllocateThroughAxom()
{
}

conduit::index_t ConduitAllocateThroughAxom::getConduitAllocatorID() const
{
  return conduitAllocatorID;
}

void *ConduitAllocateThroughAxom::internal_allocate(size_t items, size_t item_size)
{
  void *ptr = static_cast<void *>(axom::allocate<std::uint8_t>(items * item_size, axomAllocatorID));
  //std::cout << "Allocating for Conduit via axom: items=" << items << ", item_size=" << item_size << ", ptr=" << ptr << std::endl;
  return ptr;
}

void ConduitAllocateThroughAxom::internal_free(void *ptr)
{
  //std::cout << "Dellocating for Conduit via axom: ptr=" << ptr << std::endl;
  axom::deallocate(ptr);
}
#endif

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom
