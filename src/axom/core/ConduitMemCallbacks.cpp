// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/ConduitMemCallbacks.hpp"

namespace axom
{
  std::map<int, std::shared_ptr<ConduitMemCallbacks>> ConduitMemCallbacks::s_axomToInstance;
  std::map<conduit::index_t, std::shared_ptr<ConduitMemCallbacks>> ConduitMemCallbacks::s_conduitToInstance;

  const ConduitMemCallbacks& ConduitMemCallbacks::instanceForAxomId(int axomAllocId)
  {
    // This method IS NOT thread safe.

    if(s_axomToInstance.empty())
    {
      // Required one-time actions
      static auto axomMemcopy = [](void* dst, const void* src, size_t byteCount) {
        axom::copy(dst, src, byteCount);
      };
      static auto axomMemset = [](void* ptr, int value, size_t count) {
        if(axom::getAllocatorIDFromPointer(ptr) == axom::DYNAMIC_ALLOCATOR_ID)
        {
          std::memset(ptr, value, count);
        }
        else
        {
          umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
          rm.memset(ptr, value, count);
        }
      };
      conduit::utils::set_memcpy_handler(axomMemcopy);
      conduit::utils::set_memset_handler(axomMemset);
    }

    auto it = s_axomToInstance.find(axomAllocId);
    if(it == s_axomToInstance.end())
    {
      std::shared_ptr<ConduitMemCallbacks> newInstance(new ConduitMemCallbacks(axomAllocId));
      s_axomToInstance[axomAllocId] = newInstance;
      it = s_axomToInstance.insert({axomAllocId, newInstance}).first;

      auto conduitAllocId = newInstance->m_conduitId;
      assert(s_conduitToInstance.find(conduitAllocId) == s_conduitToInstance.end());
      s_conduitToInstance[conduitAllocId] = newInstance;
    }
    assert(it->first == axomAllocId);

    return *it->second;
  }

  const ConduitMemCallbacks& ConduitMemCallbacks::instanceForConduitId(conduit::index_t conduitAllocId)
  {
    // This method IS thread safe.

    auto it = s_conduitToInstance.find(conduitAllocId);
    if(it == s_conduitToInstance.end())
    {
      std::cerr << "*** Error: Axom allocator for Conduit allocator id " << conduitAllocId << " doesn't exist.  You have to register the Axom allocator first, using ConduitMemCallbacks::getAxomInstance()." << std::endl;
      axom::utilities::processAbort();
    }
    assert(it->first == conduitAllocId);

    return *it->second;
  }
}  // end namespace axom
