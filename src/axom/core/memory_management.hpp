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

#include "axom/config.hpp" // for AXOM compile-time definitions

#include <cstdlib>     // for std::malloc, std::realloc, std::free

namespace axom
{

/*!
 * \brief Allocates a chunk of memory of type T.
 * \param [in] n the number of elements to allocate.
 * \tparam T the type of pointer returned.
 * \return A pointer to the new allocation or a null pointer if allocation
 *  failed.
 */
template < typename T >
inline T* alloc( std::size_t n )
{
  return static_cast< T* >( std::malloc( n * sizeof( T ) ) );
}

/*!
 * \brief Reallocates the chunk of memory pointed to by pointer.
 * \param [in] pointer pointer to memory previously allocated with
 *  alloc or realloc, or a null pointer.
 * \param [in] n the number of elements to allocate.
 * \tparam T the type pointer points to.
 * \return A pointer to the new allocation or a null pointer if allocation
 *  failed.
 */
template < typename T >
inline T* realloc( T* pointer, std::size_t n )
{
  if ( n == 0 )
  {
    std::free( pointer );
    return nullptr;
  }

  return static_cast< T* >( std::realloc( pointer,  n * sizeof( T ) ) );
}

/*!
 * \brief Frees the chunk of memory pointed to by pointer.
 *
 * \param [in] pointer pointer to memory previously allocated with
 *  alloc or realloc or a null pointer.
 *
 *  \post pointer == nullptr
 */
template < typename T >
inline void free( T*& pointer )
{
  std::free( pointer );
  pointer = nullptr;
}

} // namespace axom

#endif /* AXOM_MEMORYMANAGEMENT_HPP_ */
