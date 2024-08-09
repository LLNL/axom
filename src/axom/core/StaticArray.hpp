// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_STATICARRAY_HPP_
#define AXOM_STATICARRAY_HPP_

#include "axom/config.hpp"  // for compile-time defines
#include "axom/core/Macros.hpp"
#include "axom/core/StackArray.hpp"

namespace axom
{
/**
 * \brief This class extends StackArray with some std::vector-like convenience methods.
 *
 * \tparam T the type of the values to hold.
 * \tparam N the number of values in the array.
 *
 */
template <typename T, int N>
class StaticArray : public StackArray<T, N>
{
public:
  AXOM_HOST_DEVICE
  constexpr axom::IndexType capacity() const
  {
    return static_cast<axom::IndexType>(N);
  }

  AXOM_HOST_DEVICE
  axom::IndexType size() const { return m_size; }

  AXOM_HOST_DEVICE
  void push_back(const T &e)
  {
    if(m_size + 1 <= capacity()) this->m_data[m_size++] = e;
  }

  AXOM_HOST_DEVICE
  void pop_back() { m_size = (m_size > 0) ? (m_size - 1) : 0; }

  AXOM_HOST_DEVICE
  void clear() { m_size = 0; }

  AXOM_HOST_DEVICE
  bool empty() const { return m_size == 0; }

  AXOM_HOST_DEVICE
  void fill(const T &e)
  {
    for(size_t i = 0; i < capacity(); i++) this->m_data[i] = e;
  }

private:
  axom::IndexType m_size {0};
};

} /* namespace axom */

#endif /* AXOM_STATICARRAY_HPP_ */
