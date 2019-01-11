/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

#ifndef AXOM_MEMORYMANAGEMENT_HPP_
#define AXOM_MEMORYMANAGEMENT_HPP_

// Axom includes
#include "axom/config.hpp" // for AXOM compile-time definitions

// Umpire includes
#ifdef AXOM_USE_UMPIRE
#include "umpire/config.hpp"
#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"
#endif

// C/C++ includes
#include <cassert>  // for assert()
#include <cstdlib>  // for std::malloc, std::realloc, std::free

namespace axom
{

/*!
 * \brief Enumerates the available memory spaces on a given system.
 *
 * \note The number of memory spaces available depends on the target system and
 *  whether Axom is compiled with CUDA and Umpire. Specifically, HOST memory is
 *  the default and is always available. If CUDA and UMPIRE support is enabled
 *  the following memory spaces are also available:
 *    * HOST_PINNED
 *    * DEVICE
 *    * DEVICE_CONSTANT
 *    * UNIFIED_MEMORY
 *
 */
enum class MemorySpace : int
{
  HOST,

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)

  HOST_PINNED,
  DEVICE,
  DEVICE_CONSTANT,
  UNIFIED_MEMORY

#endif /* defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE) */

};

/// \name Internal Data Structures
/// @{

namespace internal
{

/*!
 * \brief Holds the value to the default memory space.
 *
 * \note The default memory space is set to HOST, which is always available.
 *  The default may be changed by calling axom::setDefaultMemorySpace().
 */
static MemorySpace s_mem_space = MemorySpace::HOST;

#ifdef AXOM_USE_UMPIRE

/*!
 * \brief Maps a MemorySpace enum to the corresponding Umpire resource type.
 */
static const int umpire_type[]  =
{
  umpire::resource::Host,

#ifdef AXOM_USE_CUDA

  umpire::resource::Pinned,
  umpire::resource::Device,
  umpire::resource::Constant,
  umpire::resource::Unified

#endif /* AXOM_USE_CUDA */

};

#endif /* AXOM_USE_UMPIRE */

} /* end namspace internal */

/// @}


/// \name Memory Management Routines
/// @{

/*!
 * \brief Sets the default memory space to use. Default is set to HOST
 * \param [in] spaceId ID of the memory space to use.
 */
inline void setDefaultMemorySpace( MemorySpace spaceId );

/*!
 * \brief Returns the current default memory space used.
 * \return memSpace the current default memory space used.
 */
inline MemorySpace getDefaultMemorySpace( ) {
  return internal::s_mem_space;
};

/*!
 * \brief Allocates a chunk of memory of type T.
 *
 * \param [in] n the number of elements to allocate.
 * \param [in] spaceId the memory space where memory will be allocated
 *(optional)
 *
 * \tparam T the type of pointer returned.
 *
 * \note By default alloc() will use the current default memory space. The
 *  caller may explicitly specify the memory space to use by specifying the
 *  second, optional argument, or change the default memory space by calling
 *  axom::setDefaultMemorySpace().
 *
 * \return p pointer to the new allocation or a nullptr if allocation failed.
 *
 * \pre spaceId >= 0 && spaceId < NUM_MEMORY_SPACES
 */
template < typename T >
inline T* alloc( std::size_t n, MemorySpace spaceId=internal::s_mem_space );

/*!
 * \brief Frees the chunk of memory pointed to by the supplied pointer, p.
 * \param [in] p a pointer to memory allocated with alloc/realloc or a nullptr.
 * \post p == nullptr
 */
template < typename T >
inline void free( T*& p );

/*!
 * \brief Reallocates the chunk of memory pointed to by the supplied pointer.
 *
 * \param [in] p pointer to memory allocated with alloc/realloc, or a nullptr.
 * \param [in] n the number of elements to allocate.
 *
 * \tparam T the type pointer p points to.
 *
 * \return p pointer to the new allocation or a nullptr if allocation failed.
 */
template < typename T >
inline T* realloc( T* p, std::size_t n );

/// @}

//------------------------------------------------------------------------------
//                        IMPLEMENTATION
//------------------------------------------------------------------------------

inline void setDefaultMemorySpace( MemorySpace spaceId )
{
  internal::s_mem_space = spaceId;

#ifdef AXOM_USE_UMPIRE

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  umpire::Allocator allocator = rm.getAllocator(
    ( internal::umpire_type[ static_cast< int >( spaceId ) ] ) );

  rm.setDefaultAllocator( allocator );

#endif

}

//------------------------------------------------------------------------------
template < typename T >
inline T* alloc( std::size_t n, MemorySpace spaceId )
{
  const std::size_t numbytes = n * sizeof( T );
  T* ptr = nullptr;

#ifdef AXOM_USE_UMPIRE

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  umpire::Allocator allocator =
    rm.getAllocator( internal::umpire_type[ static_cast< int >( spaceId ) ] );

  ptr = static_cast< T* >( allocator.allocate( numbytes )  );

#else

  /* silence compiler warnings */
  static_cast< void >( spaceId );

  ptr = static_cast< T* >( std::malloc( numbytes ) );

#endif

  return ptr;
}

//------------------------------------------------------------------------------
template < typename T >
inline void free( T*& pointer )
{

#ifdef AXOM_USE_UMPIRE

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator( pointer );
  allocator.deallocate( pointer );

#else

  std::free( pointer );

#endif

  pointer = nullptr;

}

//------------------------------------------------------------------------------
template < typename T >
inline T* realloc( T* pointer, std::size_t n )
{
  if ( n==0 )
  {
    axom::free( pointer );
    return nullptr;
  }

  const std::size_t numbytes = n * sizeof( T );

#ifdef AXOM_USE_UMPIRE

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  pointer = static_cast< T* >( rm.reallocate( pointer, numbytes ) );

#else

  pointer = static_cast< T* >( std::realloc( pointer, numbytes ) );

#endif

  return pointer;
}

} // namespace axom

#endif /* AXOM_MEMORYMANAGEMENT_HPP_ */
