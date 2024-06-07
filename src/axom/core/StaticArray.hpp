// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_STATICARRAY_HPP_
#define AXOM_STATICARRAY_HPP_

#include "axom/config.hpp"           // for compile-time defines
#include "axom/core/Macros.hpp"      // for axom macros
#include "axom/core/StackArray.hpp"  // for StackArray
#include "axom/core/Types.hpp"       // for axom types

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
 */
template <typename T, int N>
struct StaticArray : public StackArray<T, N>
{
  AXOM_HOST_DEVICE int size() const { return m_size; }
  AXOM_HOST_DEVICE void push_back(const T& obj)
  {
    assert(m_size < N);
    if(m_size < N)
    {
      StackArray<T, N>::m_data[m_size++] = obj;
    }
  }

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