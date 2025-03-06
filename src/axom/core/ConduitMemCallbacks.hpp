// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file ConduitMemCallbacks.hpp
 *
 * \brief   Call-backs for using Axom memory mamagement in Conduit.
 *
 ******************************************************************************
 */

#ifndef SIDRE_CONDUITMEMCALLBACKS_HPP_
#define SIDRE_CONDUITMEMCALLBACKS_HPP_

// Standard C++ headers
#include <string>
#include <set>

#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"
#include "conduit_node.hpp"
#include "conduit_utils.hpp"

namespace axom
{

/*!
  @brief Object to handle Conduit memory operations by delegating to Axom.

  This class has no public constructor.  Use getInstance(int axomAllocId)
  to access the instance for a specific Axom allocator id.  The construction
  registers the appropriate callbacks with Conduit, including the required
  memset and memcopy callbacks.

  Examples for setting Conduit allocator ids when you have Axom allocator ids.
  @code{.cpp}
    void foo(conduit::Node& n, int axomAllocId) {
      n.set_allocator(axomAllocIdToConduit(axomAllocId));
    }

    void bar(conduit::Node& n, int axomAllocId) {
      const auto& instance = getInstance(axomAllocId);
      assert(instance.axomId() == axomAllocId);
      n.set_allocator(instance.conduitId());
    }
  @endcode
*/
struct ConduitMemCallbacks {
  //!@brief Return the Axom allocator id.
  int axomId() const { return m_axomId; }

  //!@brief Return the Conduit allocator id coresponding to axomId().
  conduit::index_t conduitId() const { return m_conduitId; }

  /*!
    @brief Convert an Axom allocator id to Conduit, registering
    a new Conduit allocator if needed.
  */
  static conduit::index_t axomAllocIdToConduit(int axomAllocId)
  {
    return getInstance(axomAllocId).conduitId();
  }

  //!@brief Return the instance for the given Axom allocator id.
  static const ConduitMemCallbacks& getInstance(int axomAllocId)
  {
    // This method is not thread safe.

    // Mapping from Axom allocator to an instance.
    static std::map<int, std::shared_ptr<ConduitMemCallbacks>> s_handlers;

    if(s_handlers.empty())
    {
      // Required one-time actions
      static auto axomMemcopy = [](void* dst, const void* src, size_t byteCount) {
                                  axom::copy(dst, src, byteCount);
                                };
      static auto axomMemset = [](void* ptr, int value, size_t count) {
                                 if (axom::getAllocatorIDFromPointer(ptr)
                                     == axom::DYNAMIC_ALLOCATOR_ID) {
                                   std::memset(ptr, value, count);
                                 } else {
                                   umpire::ResourceManager& rm =
                                     umpire::ResourceManager::getInstance();
                                   rm.memset(ptr, value, count);
                                 }
                               };
      conduit::utils::set_memcpy_handler(axomMemcopy);
      conduit::utils::set_memset_handler(axomMemset);
    }

    auto it = s_handlers.find(axomAllocId);
    if(it == s_handlers.end())
    {
      it =
        s_handlers.emplace(axomAllocId, new ConduitMemCallbacks(axomAllocId)).first;
    }
    assert(it->first == axomAllocId);

    return *it->second;
  }

  ~ConduitMemCallbacks() { }

private:
  //!@brief Axom's allocator id.
  int m_axomId;
  //!@brief Conduit's allocator id equivalent to m_axomId.
  conduit::index_t m_conduitId;

#if 0
  // Conduit will support std::function in the near future.
  using AllocatorCallback = std::function<void*(size_t, size_t)>;
  using DeallocCallback = std::function<void(void*)>;
#else
  typedef void* (AllocatorCallback)(size_t, size_t);
  typedef void (DeallocCallback)(void *);
#endif

  AllocatorCallback *m_allocCallback;
  DeallocCallback *m_deallocCallback;

  ConduitMemCallbacks() = delete;

  /*!
    @brief Constructor creates allocator/deallocator function and registers
    them with Conduit.
  */
  ConduitMemCallbacks(int axomAllocId)
    : m_axomId(axomAllocId)
  {
    using conduit::utils::register_allocator;
    auto deallocator = [](void* ptr) {
      char* cPtr = (char*)(ptr);
      axom::deallocate<char>(cPtr);
    };
    m_deallocCallback = deallocator;
    /*
      Note: Once Conduit allows the callbacks as std::function types,
      we can use a single allocator, eliminating the need for these
      if-else blocks.
    */
    if(axomAllocId == 0)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 0);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 1)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 1);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 2)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 2);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 3)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 3);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 4)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 4);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else if(axomAllocId == 5)
    {
      m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
        void* ptr = axom::allocate<char>(itemCount * itemByteSize, 5);
        return ptr;
      };
      m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
    }
    else
    {
      std::cerr<<
        "*** Work-around for conduit::utils::register_allocator needs case for "
        "axomAllocId = " << std::to_string(axomAllocId) <<
        ".  Please add it to ConduitMemCallbacks.hpp.";
    }
  }

};

} /* end namespace axom */

#endif /* SIDRE_CONDUITMEMCALLBACKS_HPP_ */
