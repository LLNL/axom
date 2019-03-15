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

/// \name Memory Management Routines
/// @{

#ifdef AXOM_USE_UMPIRE

/*!
 * \brief Sets the default memory space to use. Default is set to HOST
 * \param [in] spaceId ID of the memory space to use.
 */
inline void setDefaultAllocator( umpire::Allocator allocator )
{
  umpire::ResourceManager::getInstance().setDefaultAllocator( allocator );
}

/*!
 * \brief Returns the current default memory space used.
 */
inline umpire::Allocator getDefaultAllocator()
{
  return umpire::ResourceManager::getInstance().getDefaultAllocator();
}

/*!
 * \brief Returns the umpire allocator associated with the given ID.
 * \param [in] allocatorID the ID of the allocator to get.
 */
inline umpire::Allocator getAllocator( int allocatorID )
{
  return umpire::ResourceManager::getInstance().getAllocator( allocatorID );
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
 * \note By default alloc() will use the current default memory space. The
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
inline T* alloc( std::size_t n, umpire::Allocator allocator=getDefaultAllocator() ) noexcept;
#else
inline T* alloc( std::size_t n ) noexcept;
#endif

/*!
 * \brief Frees the chunk of memory pointed to by the supplied pointer, p.
 * \param [in] p a pointer to memory allocated with alloc/realloc or a nullptr.
 * \post p == nullptr
 */
template < typename T >
inline void free( T*& p ) noexcept;

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
inline T* realloc( T* p, std::size_t n ) noexcept;

/// @}

//------------------------------------------------------------------------------
//                        IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#ifdef AXOM_USE_UMPIRE

template < typename T >
inline T* alloc( std::size_t n, umpire::Allocator allocator ) noexcept
{
  if ( n == 0 ) return nullptr;

  const std::size_t numbytes = n * sizeof( T );
  return static_cast< T* >( allocator.allocate( numbytes )  );
}

#else

template < typename T >
inline T* alloc( std::size_t n ) noexcept
{
  if ( n == 0 ) return nullptr;

  const std::size_t numbytes = n * sizeof( T );
  return static_cast< T* >( std::malloc( numbytes )  );
}

#endif

//------------------------------------------------------------------------------
template < typename T >
inline void free( T*& pointer ) noexcept
{
  if ( pointer == nullptr ) return;

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
inline T* realloc( T* pointer, std::size_t n ) noexcept
{
  if ( n == 0 )
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
