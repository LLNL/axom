// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mir/utilities.hpp"

namespace axom
{
namespace mir
{
namespace utilities
{

AXOM_HOST_DEVICE
std::uint64_t hash_bytes(const std::uint8_t *data, std::uint32_t length)
{
  std::uint32_t hash = 0;

  // Build the length into the hash.
  const auto ldata = reinterpret_cast<const std::uint8_t *>(&length);
  for(int e = 0; e < 4; e++)
  {
    hash += ldata[e];
    hash += hash << 10;
    hash ^= hash >> 6;
  }

  std::uint32_t hashr = hash;
  for(std::uint32_t i = 0; i < length; i++)
  {
    hash += data[i];
    hash += hash << 10;
    hash ^= hash >> 6;

    hashr += data[length - 1 - i];
    hashr += hashr << 10;
    hashr ^= hashr >> 6;
  }
  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;

  hashr += hashr << 3;
  hashr ^= hashr >> 11;
  hashr += hashr << 15;

  return (static_cast<std::uint64_t>(hash) << 32) | hashr;
}

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
  std::cout << "Allocating for Conduit via axom: items=" << items << ", item_size=" << item_size << ", ptr=" << ptr << std::endl;
  return ptr;
}

void ConduitAllocateThroughAxom::internal_free(void *ptr)
{
  std::cout << "Dellocating for Conduit via axom: ptr=" << ptr << std::endl;
  axom::deallocate(ptr);
}

} // end namespace utilities
} // end namespace mir
} // end namespace axom
