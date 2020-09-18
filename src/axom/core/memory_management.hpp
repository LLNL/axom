// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MEMORYMANAGEMENT_HPP_
#define AXOM_MEMORYMANAGEMENT_HPP_

// Axom includes
#include "axom/config.hpp"  // for AXOM compile-time definitions

// Umpire includes
#ifdef AXOM_USE_UMPIRE
  #include "umpire/config.hpp"
  #include "umpire/ResourceManager.hpp"
  #include "umpire/op/MemoryOperationRegistry.hpp"
#else
  #include <cstring>  // for std::memcpy
  #include <cstdlib>  // for std::malloc, std::realloc, std::free
#endif

namespace axom
{
#ifndef AXOM_USE_UMPIRE
  constexpr int DEFAULT_ALLOCATOR_ID = 0;
#endif

constexpr int INVALID_ALLOCATOR_ID = -1;

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

#endif

/*!
 * \brief Sets the default memory space to use. Default is set to HOST
 * \param [in] allocatorID ID of the allocator to use.
 */
inline void setDefaultAllocator(int allocatorID)
{
#ifdef AXOM_USE_UMPIRE
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator(allocatorID);
  rm.setDefaultAllocator(allocator);
#else
  static_cast<void>(allocatorID);  // silence compiler warnings
#endif
}

/*!
 * \brief Returns the current default memory space used.
 * \note If Umpire is used, the corresponding umpire allocator can be retrieved
 * by:
 *  <code>
 *    umpire::Allocator alloc =
 * umpire::ResourceManager::getInstance().getAllocator( allocID );
 *  </code>
 */
inline int getDefaultAllocatorID()
{
#ifdef AXOM_USE_UMPIRE
  return umpire::ResourceManager::getInstance().getDefaultAllocator().getId();
#else
  return axom::DEFAULT_ALLOCATOR_ID;
#endif
}

/*!
 * \brief Allocates a chunk of memory of type T.
 *
 * \param [in] n the number of elements to allocate.
 * \param [in] allocator the Umpire allocator to use
 *(optional)
 *
 * \tparam T the type of pointer returned.
 *
 * \note By default allocate() will use the current default memory space. The
 *  caller may explicitly specify the memory space to use by specifying the
 *  second, optional argument, or change the default memory space by calling
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
 *
 * \tparam T the type pointer p points to.
 *
 * \return p pointer to the new allocation or a nullptr if allocation failed.
 *
 * \note When n == 0, this function returns a valid pointer (of size 0) in the
 * current allocator's memory space. This follows the semantics of
 * Umpire's reallocate function.
 */
template <typename T>
inline T* reallocate(T* p, std::size_t n) noexcept;

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
inline void copy(void* dst, void* src, std::size_t numbytes) noexcept;

/// @}

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
  static_cast<void>(allocID);  // silence compiler warnings
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
inline T* reallocate(T* pointer, std::size_t n) noexcept
{
  const std::size_t numbytes = n * sizeof(T);

#if defined(AXOM_USE_UMPIRE) && !defined(UMPIRE_VERSION_MAJOR)

  // Workaround for bug in Umpire's handling on reallocate(0)
  // Fixed in Umpire PR #292 (after v1.1.0)

  // NOTE: The UMPIRE_VERSION_MAJOR macro was added in umpire-v2.0.0. If the
  // macro is not defined, we assume that the Umpire version is less than 2.0.0
  // and that the workaround is needed.
  if(n == 0)
  {
    axom::deallocate<T>(pointer);
    pointer = axom::allocate<T>(0);
    return pointer;
  }

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  // Workaround for bug in Umpire's handling of reallocate
  // called on a zero-sized allocation
  // Fixed in Umpire PR #292 (after v1.1.0)
  if(pointer != nullptr)
  {
    auto* allocRecord = rm.findAllocationRecord(pointer);
    if(allocRecord && allocRecord->size == 0)
    {
      axom::deallocate<T>(pointer);
      pointer = axom::allocate<T>(n);
      return pointer;
    }
  }

  pointer = static_cast<T*>(rm.reallocate(pointer, numbytes));

#elif defined(AXOM_USE_UMPIRE) && \
      ( (UMPIRE_VERSION_MAJOR == 2) && (UMPIRE_VERSION_MINOR >= 1) ) || \
      (UMPIRE_VERSION_MAJOR > 2 )

  // Umpire 2.1.0 and above handles reallocate(0) natively
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  pointer = static_cast<T*>(rm.reallocate(pointer, numbytes));

#else

  pointer = static_cast<T*>(std::realloc(pointer, numbytes));

  // Consistently handle realloc(0) for std::realloc to match Umpire's behavior
  if(n == 0 && pointer == nullptr)
  {
    pointer = axom::allocate<T>(0);
  }

#endif

  return pointer;
}

inline void copy(void* dst, void* src, std::size_t numbytes) noexcept
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
    dstRecord   = const_cast<AllocationRecord*>(rm.findAllocationRecord(dst));
    dstStrategy = dstRecord->strategy;
  }

  if(rm.hasAllocator(src))
  {
    srcRecord   = const_cast<AllocationRecord*>(rm.findAllocationRecord(src));
    srcStrategy = srcRecord->strategy;
  }

  auto op = op_registry.find("COPY", srcStrategy, dstStrategy);
  op->transform(src, &dst, srcRecord, dstRecord, numbytes);
#else
  std::memcpy(dst, src, numbytes);
#endif
}

}  // namespace axom

#endif /* AXOM_MEMORYMANAGEMENT_HPP_ */
