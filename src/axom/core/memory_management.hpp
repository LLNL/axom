// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MEMORYMANAGEMENT_HPP_
#define AXOM_MEMORYMANAGEMENT_HPP_

// Axom includes
#include "axom/config.hpp" // for AXOM compile-time definitions

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

constexpr int INVALID_ALLOCATOR_ID = -1;

/// \name Memory Management Routines
/// @{

#ifdef AXOM_USE_UMPIRE

/*!
 * \brief Returns the ID of the predefined allocator for a given resource.
 * \param [in] resource_type the Umpire resource type
 * \return ID the id of the predefined umpire allocator.
 */
inline int getResourceAllocatorID(
  umpire::resource::MemoryResourceType resource_type )
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator alloc     = rm.getAllocator( resource_type );
  return alloc.getId();
}

/*!
 * \brief Returns the umpire allocator associated with the given ID.
 * \param [in] allocatorID the ID of the allocator to get.
 */
inline umpire::Allocator getAllocator( int allocatorID )
{
  return umpire::ResourceManager::getInstance().getAllocator( allocatorID );
}

/*!
 * \brief Sets the default memory space to use. Default is set to HOST
 * \param [in] allocator the umpire::Allocator to make default.
 */
inline void setDefaultAllocator( umpire::Allocator allocator )
{
  umpire::ResourceManager::getInstance().setDefaultAllocator( allocator );
}

/*!
 * \brief Sets the default memory space to use. Default is set to HOST
 * \param [in] allocatorID ID of the umpire::Allocator to use.
 */
inline void setDefaultAllocator( int allocatorID )
{
  setDefaultAllocator( getAllocator( allocatorID ) );
}

/*!
 * \brief Returns the current default memory space used.
 */
inline umpire::Allocator getDefaultAllocator()
{
  return umpire::ResourceManager::getInstance().getDefaultAllocator();
}

#endif

/*!
 * \brief Allocates a chunk of memory of type T.
 *
 * \param [in] n the number of elements to allocate.
 * \param [in] spaceId the memory space where memory will be allocated
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
 *
 * \pre spaceId >= 0 && spaceId < NUM_MEMORY_SPACES
 */
template < typename T >
#ifdef AXOM_USE_UMPIRE
inline T* allocate( std::size_t n,
                    umpire::Allocator allocator=
                      getDefaultAllocator() ) noexcept;
#else
inline T* allocate( std::size_t n ) noexcept;
#endif

/*!
 * \brief Frees the chunk of memory pointed to by the supplied pointer, p.
 * \param [in/out] p a pointer to memory allocated with allocate/reallocate or a
 * nullptr.
 * \post p == nullptr
 */
template < typename T >
inline void deallocate( T*& p ) noexcept;

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
 */
template < typename T >
inline T* reallocate( T* p, std::size_t n ) noexcept;

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
inline void copy( void* dst, void* src, std::size_t numbytes ) noexcept;

/// @}

//------------------------------------------------------------------------------
//                        IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#ifdef AXOM_USE_UMPIRE

template < typename T >
inline T* allocate( std::size_t n, umpire::Allocator allocator ) noexcept
{
  const std::size_t numbytes = n * sizeof( T );
  return static_cast< T* >( allocator.allocate( numbytes )  );
}

#else

template < typename T >
inline T* allocate( std::size_t n ) noexcept
{
  const std::size_t numbytes = n * sizeof( T );
  return static_cast< T* >( std::malloc( numbytes )  );
}

#endif

//------------------------------------------------------------------------------
template < typename T >
inline void deallocate( T*& pointer ) noexcept
{
  if ( pointer == nullptr )
    return;

#ifdef AXOM_USE_UMPIRE

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  rm.deallocate( pointer );

#else

  std::free( pointer );

#endif

  pointer = nullptr;
}

//------------------------------------------------------------------------------
template < typename T >
inline T* reallocate( T* pointer, std::size_t n ) noexcept
{
  const std::size_t numbytes = n * sizeof( T );

#ifdef AXOM_USE_UMPIRE

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  pointer = static_cast< T* >( rm.reallocate( pointer, numbytes ) );

#else

  pointer = static_cast< T* >( std::realloc( pointer, numbytes ) );

#endif

  return pointer;
}

inline void copy( void* dst, void* src, std::size_t numbytes ) noexcept
{
#ifdef AXOM_USE_UMPIRE
  umpire::ResourceManager & rm = umpire::ResourceManager::getInstance();
  umpire::op::MemoryOperationRegistry & op_registry =
    umpire::op::MemoryOperationRegistry::getInstance();

  auto dstStrategy = rm.getAllocator( "HOST" ).getAllocationStrategy();
  auto srcStrategy = dstStrategy;

  if (rm.hasAllocator(dst))
  {
    dstStrategy = rm.findAllocationRecord( dst )->strategy;
  }

  if (rm.hasAllocator(src))
  {
    srcStrategy = rm.findAllocationRecord( src )->strategy;
  }

  auto op = op_registry.find( "COPY", srcStrategy, dstStrategy );
  op->transform( src, &dst, nullptr, nullptr, numbytes );
#else
  std::memcpy( dst, src, numbytes );
#endif
}

} // namespace axom

#endif /* AXOM_MEMORYMANAGEMENT_HPP_ */
