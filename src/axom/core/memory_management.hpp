// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MEMORYMANAGEMENT_HPP_
#define AXOM_MEMORYMANAGEMENT_HPP_

// Axom includes
#include "axom/config.hpp"  // for AXOM compile-time definitions
#include "axom/core/Macros.hpp"

// Umpire includes
#ifdef AXOM_USE_UMPIRE
  #include "umpire/config.hpp"
  #include "umpire/ResourceManager.hpp"
  #include "umpire/op/MemoryOperationRegistry.hpp"
  #include "umpire/resource/MemoryResourceTypes.hpp"
  #include "umpire/strategy/QuickPool.hpp"
#else
  #include <cstring>  // for std::memcpy
  #include <cstdlib>  // for std::malloc, std::realloc, std::free
#endif

namespace axom
{
constexpr int INVALID_ALLOCATOR_ID = -1;

// _memory_space_start
/*! 
 * \brief Memory spaces supported by Array-like types
 *
 * This abstraction is not implemented using Umpire's MemoryResourceType enum
 * in order to also include a "Dynamic" option as a default template parameter
 * for Array-like types
 */
enum class MemorySpace
{
  Dynamic,
#ifdef AXOM_USE_UMPIRE
  Host,
  Device,
  Unified,
  Pinned,
  Constant
#endif
};
// _memory_space_end

// _memory_management_routines_start
/// \name Memory Management Routines
/// @{

#ifdef AXOM_USE_UMPIRE

/*!
 * \brief Returns the ID of the predefined allocator for a given resource.
 * \param [in] resource_type the Umpire resource type
 * \return ID the id of the predefined umpire allocator.
 */
inline int getUmpireResourceAllocatorID(
  umpire::resource::MemoryResourceType resource_type)
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator alloc = rm.getAllocator(resource_type);
  return alloc.getId();
}

/*!
 * \brief Sets the default memory allocator to use.
 * \param [in] resource_type the Umpire resource type
 */
inline void setDefaultAllocator(umpire::resource::MemoryResourceType resource_type)
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator(resource_type);
  rm.setDefaultAllocator(allocator);
}

#endif

/*!
 * \brief Sets the default memory allocator to use.
 * \param [in] allocId the Umpire allocator id
 * 
 * \note This function has no effect when Axom is not compiled with Umpire.
 */
inline void setDefaultAllocator(int allocId)
{
#ifdef AXOM_USE_UMPIRE
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator(allocId);
  rm.setDefaultAllocator(allocator);
#else
  AXOM_UNUSED_VAR(allocId);
#endif
}

/*!
 * \brief Returns the ID of the current default allocator.
 * \return ID the ID of the current default allocator.
 * \post ID != INVALID_ALLOCATOR_ID
 */
inline int getDefaultAllocatorID()
{
#ifdef AXOM_USE_UMPIRE
  return umpire::ResourceManager::getInstance().getDefaultAllocator().getId();
#else
  return 0;
#endif
}

/*!
 * \brief Allocates a chunk of memory of type T.
 *
 * \param [in] n the number of elements to allocate.
 * \param [in] allocID the Umpire allocator to use (optional)
 *
 * \tparam T the type of pointer returned.
 *
 * \note By default allocate() will use the current default allocator. The
 *  caller may explicitly specify a different allocator to use by supplying the
 *  second, optional argument, or change the default allocator by calling
 *  axom::setDefaultAllocator().
 *
 * \return p pointer to the new allocation or a nullptr if allocation failed.
 */
template <typename T>
inline T* allocate(std::size_t n, int allocID = getDefaultAllocatorID()) noexcept;

/*!
 * \brief Frees the chunk of memory pointed to by the supplied pointer, p.
 * \param [in/out] p a pointer to memory allocated with allocate/reallocate or a
 * nullptr.
 * \post p == nullptr
 */
template <typename T>
inline void deallocate(T*& p) noexcept;

/*!
 * \brief Reallocates the chunk of memory pointed to by the supplied pointer.
 *
 * \param [in] p pointer to memory allocated with allocate/reallocate, or a
 * nullptr.
 * \param [in] n the number of elements to allocate.
 * \param [in] allocID the ID of the allocator to use if pointer is null
 * (optional)
 *
 * \tparam T the type pointer p points to.
 *
 * \return p pointer to the new allocation or a nullptr if allocation failed.
 *
 * \note When n == 0, this function returns a valid pointer (of size 0) in the
 * current allocator's memory space. This follows the semantics of
 * Umpire's reallocate function.
 * \note When p is a null pointer, allocID is used to allocate the data.
 * Otherwise, it is unused.
 */
template <typename T>
inline T* reallocate(T* p,
                     std::size_t n,
                     int allocID = getDefaultAllocatorID()) noexcept;

/*!
 * \brief Copies memory from the source to the destination.
 *
 * \param [in/out] dst the destination to copy to.
 * \param [in] src the source to copy from.
 * \param [in] numbytes the number of bytes to copy.
 *
 * \note When using Umpire if either src or dst is not registered with the
 *  ResourceManager then the default host allocation strategy is assumed for
 *  that pointer.
 */
inline void copy(void* dst, const void* src, std::size_t numbytes) noexcept;

/// @}
// _memory_management_routines_end

//------------------------------------------------------------------------------
//                        IMPLEMENTATION
//------------------------------------------------------------------------------

template <typename T>
inline T* allocate(std::size_t n, int allocID) noexcept
{
  const std::size_t numbytes = n * sizeof(T);

#ifdef AXOM_USE_UMPIRE

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator(allocID);
  return static_cast<T*>(allocator.allocate(numbytes));

#else
  AXOM_UNUSED_VAR(allocID);
  return static_cast<T*>(std::malloc(numbytes));
#endif
}
//------------------------------------------------------------------------------
template <typename T>
inline void deallocate(T*& pointer) noexcept
{
  if(pointer == nullptr) return;

#ifdef AXOM_USE_UMPIRE

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  rm.deallocate(pointer);

#else

  std::free(pointer);

#endif

  pointer = nullptr;
}

//------------------------------------------------------------------------------
template <typename T>
inline T* reallocate(T* pointer, std::size_t n, int allocID) noexcept
{
  const std::size_t numbytes = n * sizeof(T);

#if defined(AXOM_USE_UMPIRE)

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  if(pointer == nullptr)
  {
    pointer = axom::allocate<T>(n, allocID);
  }
  else
  {
    pointer = static_cast<T*>(rm.reallocate(pointer, numbytes));
  }

#else

  pointer = static_cast<T*>(std::realloc(pointer, numbytes));

  // Consistently handle realloc(0) for std::realloc to match Umpire's behavior
  if(n == 0 && pointer == nullptr)
  {
    pointer = axom::allocate<T>(0);
  }

  AXOM_UNUSED_VAR(allocID);
#endif

  return pointer;
}

//------------------------------------------------------------------------------
inline void copy(void* dst, const void* src, std::size_t numbytes) noexcept
{
#ifdef AXOM_USE_UMPIRE
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::op::MemoryOperationRegistry& op_registry =
    umpire::op::MemoryOperationRegistry::getInstance();

  auto dstStrategy = rm.getAllocator("HOST").getAllocationStrategy();
  auto srcStrategy = dstStrategy;

  using AllocationRecord = umpire::util::AllocationRecord;
  AllocationRecord* dstRecord = nullptr;
  AllocationRecord* srcRecord = nullptr;

  if(rm.hasAllocator(dst))
  {
    dstRecord = const_cast<AllocationRecord*>(rm.findAllocationRecord(dst));
    dstStrategy = dstRecord->strategy;
  }

  if(rm.hasAllocator(const_cast<void*>(src)))
  {
    srcRecord = const_cast<AllocationRecord*>(
      rm.findAllocationRecord(const_cast<void*>(src)));
    srcStrategy = srcRecord->strategy;
  }

  auto op = op_registry.find("COPY", srcStrategy, dstStrategy);
  op->transform(const_cast<void*>(src), &dst, srcRecord, dstRecord, numbytes);
#else
  std::memcpy(dst, src, numbytes);
#endif
}

namespace detail
{
/// \brief Translates between the MemorySpace enum and Umpire allocator IDs
template <MemorySpace SPACE>
inline int getAllocatorID();

template <>
inline int getAllocatorID<MemorySpace::Dynamic>()
{
  return axom::getDefaultAllocatorID();
}

inline MemorySpace getAllocatorSpace(int allocatorId)
{
#ifdef AXOM_USE_UMPIRE
  using ump_res_type = typename umpire::MemoryResourceTraits::resource_type;

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  auto umpResType =
    rm.getAllocator(allocatorId).getAllocationStrategy()->getTraits().resource;
  switch(umpResType)
  {
  case ump_res_type::host:
    return MemorySpace::Host;
  case ump_res_type::device:
    return MemorySpace::Device;
  case ump_res_type::device_const:
    return MemorySpace::Constant;
  case ump_res_type::pinned:
    return MemorySpace::Pinned;
  case ump_res_type::um:
    return MemorySpace::Unified;
  default:
    return MemorySpace::Dynamic;
  }
#else
  AXOM_UNUSED_VAR(allocatorId);
  return MemorySpace::Dynamic;
#endif
}

#ifdef AXOM_USE_UMPIRE

template <>
inline int getAllocatorID<MemorySpace::Host>()
{
  return axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Host);
}

template <>
inline int getAllocatorID<MemorySpace::Device>()
{
  return axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Device);
}

template <>
inline int getAllocatorID<MemorySpace::Unified>()
{
  return axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);
}

template <>
inline int getAllocatorID<MemorySpace::Pinned>()
{
  return axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Pinned);
}

template <>
inline int getAllocatorID<MemorySpace::Constant>()
{
  return axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Constant);
}

#endif

}  // namespace detail

}  // namespace axom

#endif /* AXOM_MEMORYMANAGEMENT_HPP_ */
