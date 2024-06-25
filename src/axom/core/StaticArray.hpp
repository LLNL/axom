// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_STATICARRAY_HPP_
#define AXOM_STATICARRAY_HPP_

#include "axom/config.hpp"           // for compile-time defines
#include "axom/core/Macros.hpp"      // for axom macros
#include "axom/core/StackArray.hpp"  // for StackArray

namespace axom
{
/*!
 * \accelerated
 * \class StaticArray
 *
 * \brief Provides a wrapper for a StackArray with an additional size member.
 *
 * \tparam T the type of the values to hold.
 * \tparam N the number of values in the array.
 *
 * \note Type \a T must be default-constructible on device for device
 *       execution.
 */
template <typename T, int N>
struct StaticArray : public StackArray<T, N>
{
  /*!
   * \brief Returns the size of the static array
   *
   * \return The size of the static array
   */
  AXOM_HOST_DEVICE int size() const { return m_size; }

  /*!
   * \brief Pushes an object to the back of the static array
   *
   * \param [in] obj the object to be added to the back.
   *
   * \note The number of push_backs must not exceed N,
   *       the max number of values in the array.
   *
   * \note If the static array is full, push_back
   *       will not modify the static array.
   */
  AXOM_HOST_DEVICE void push_back(const T& obj)
  {
    assert(m_size < N);
    if(m_size < N)
    {
      StackArray<T, N>::m_data[m_size++] = obj;
    }
  }

  /*!
   * \brief Clears the data from the static array
   */
  AXOM_HOST_DEVICE void clear()
  {
    for(T& datum : StackArray<T, N>::m_data)
    {
      datum = T();
    }
    m_size = 0;
  }

private:
  int m_size {0};
};

} /* namespace axom */

#endif /* AXOM_STATICARRAY_HPP_ */
