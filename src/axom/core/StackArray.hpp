// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_STACKARRAY_HPP_
#define AXOM_STACKARRAY_HPP_

#include "axom/config.hpp"       // for compile-time defines
#include "axom/core/Macros.hpp"  // for axom macros
#include "axom/core/Types.hpp"   // for axom types

namespace axom
{
/*!
 * \class StackArray
 *
 * \brief Provides a wrapper for a compile time sized array, similar to
 *  std::array. This class is needed because NVCC doesn't caputure standard
 *  stack arrays in device lambdas. Furthermore we can't use std::array becuase
 *  it is not host-device decorated.
 *
 * \tparam T the type of the values to hold.
 * \tparam N the number of values in the array.
 */
template <typename T, int N>
struct StackArray
{
  /*!
   * \brief Accessor, returns a reference to the value at the given index.
   *
   * \param [in] i the index to access.
   */
  /// @{

  AXOM_HOST_DEVICE
  T& operator[](IndexType i) noexcept { return m_data[i]; }

  AXOM_HOST_DEVICE
  constexpr const T& operator[](IndexType i) const noexcept
  {
    return m_data[i];
  }

  /// @}

  /*!
   * \brief User defined conversion to a raw pointer.
   */
  /// @{

  AXOM_HOST_DEVICE operator T*() noexcept { return &m_data[0]; }

  AXOM_HOST_DEVICE
  constexpr operator const T*() const noexcept { return &m_data[0]; }

  /// @}

  T m_data[N];
};

} /* namespace axom */

#endif /* AXOM_STACKARRAY_HPP_ */
