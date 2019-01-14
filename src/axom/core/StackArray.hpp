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

#ifndef AXOM_STACKARRAY_HPP_
#define AXOM_STACKARRAY_HPP_

#include "axom/config.hpp"                   // for compile-time defines
#include "axom/core/Macros.hpp"              // for axom macros

namespace axom
{

template< typename T, int N >
struct StackArray
{
  /*!
   * \brief Accessor, returns a reference to the value at the given index.
   *
   * \param [in] i the index to access.
   */
  /// @{

  AXOM_HOST_DEVICE
  T& operator[](int i) noexcept
  { return m_data[ i ]; }

  AXOM_HOST_DEVICE 
  constexpr const T& operator[](int i) const noexcept
  { return m_data[i]; }

  /// @}

  /*!
   * \brief User defined conversion to a raw pointer.
   */
  /// @{

  AXOM_HOST_DEVICE operator T* () noexcept
  { return &m_data[0]; }

  AXOM_HOST_DEVICE 
  constexpr operator const T* () const noexcept
  { return &m_data[0]; }

  /// @}

  T m_data[  N  ];
};

/*!
 * \brief Creates a StackArray that copies the data pointed to by source.
 *  This is a free method so that StackArray doesn't define a constructor.
 *
 * \param [in] source the data to copy.
 */
template< typename T, int N >
StackArray< T, N > createStackArray( const T * source )
{
  StackArray< T, N > arr;
  for ( int i = 0; i < N; ++i )
  {
    arr[ i ] = source[ i ];
  }

  return arr;
}

} /* namespace axom */

#endif /* AXOM_STACKARRAY_HPP_ */
