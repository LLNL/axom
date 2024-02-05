// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_ARRAYINDEXER_HPP_
#define QUEST_ARRAYINDEXER_HPP_

#include "axom/slic.hpp"
#include "axom/core/StackArray.hpp"
#include "axom/core/numerics/matvecops.hpp"

namespace axom
{
/*!
  @brief Indexing into a multidimensional structured array.

  Supports row-major and column-major ordering and arbitrary
  permutations of the ordering.
*/
template <typename T, int DIM>
class ArrayIndexer
{
  axom::StackArray<T, DIM> m_strides;
  axom::StackArray<std::uint16_t, DIM> m_slowestDirs;

public:
  /*!
    @brief Constructor for row- or column major indexing.
    @param [in] lengths Lengths of the array
    @param [in] order: c is column major; r is row major.
  */
  ArrayIndexer(const axom::StackArray<T, DIM>& lengths, char order)
  {
    initialize(lengths, order);
  }

  //!@brief Constructor for arbitrary-stride indexing.
  ArrayIndexer(const axom::StackArray<T, DIM>& strides) : m_strides(strides)
  {
    initialize(strides);
  }

  /*!
    @brief Initialize for row- or column major indexing.
    @param [in] shape Shape of the array
    @param [in] order: c is column major; r is row major.
  */
  inline AXOM_HOST_DEVICE void initialize(const axom::StackArray<T, DIM>& shape,
                                          char order)
  {
    SLIC_ASSERT(order == 'c' || order == 'r');
    if(order == 'r')
    {
      for(int d = 0; d < DIM; ++d)
      {
        m_slowestDirs[d] = DIM - 1 - d;
      }
      m_strides[0] = 1;
      for(int d = 1; d < DIM; ++d)
      {
        m_strides[d] = m_strides[d - 1] * shape[d - 1];
      }
    }
    else
    {
      for(int d = 0; d < DIM; ++d)
      {
        m_slowestDirs[d] = d;
      }
      m_strides[DIM - 1] = 1;
      for(int d = DIM - 2; d >= 0; --d)
      {
        m_strides[d] = m_strides[d + 1] * shape[d + 1];
      }
    }
    SLIC_ASSERT((DIM == 1 && getOrder() == ('r' | 's')) || (getOrder() == order));
  }

  //!@brief Initialize for arbitrary-stride indexing.
  inline AXOM_HOST_DEVICE void initialize(const axom::StackArray<T, DIM>& strides)
  {
    m_strides = strides;
    for(int d = 0; d < DIM; ++d)
    {
      m_slowestDirs[d] = d;
    }
    for(int s = 0; s < DIM; ++s)
    {
      for(int d = s; d < DIM; ++d)
      {
        if(m_strides[m_slowestDirs[s]] < m_strides[m_slowestDirs[d]])
        {
          // Swap values.
          auto tmp = m_slowestDirs[s];
          m_slowestDirs[s] = m_slowestDirs[d];
          m_slowestDirs[d] = tmp;
        }
      }
    }
  }

  //!@brief Index directions, ordered from slowest to fastest.
  inline AXOM_HOST_DEVICE const axom::StackArray<std::uint16_t, DIM>& slowestDirs() const
  {
    return m_slowestDirs;
  }

  //!@brief Strides.
  inline AXOM_HOST_DEVICE const axom::StackArray<axom::IndexType, DIM>& strides() const
  {
    return m_strides;
  }

  /*!
    @brief Get the stride order (row- or column-major, or something else).

    @return 'r' or 'c' for row- or column major, '\0' for neither,
    or if DIM == 1, the value of 'r' | 'c'.
  */
  inline AXOM_HOST_DEVICE char getOrder() const
  {
    char order = 'r' | 'c';
    for(int d = 0; d < DIM - 1; ++d)
    {
      order &= m_slowestDirs[d] < m_slowestDirs[d + 1] ? 'c' : 'r';
    }
    return order;
  }

  //!@brief Convert multidimensional index to flat index.
  inline AXOM_HOST_DEVICE T toFlatIndex(const axom::StackArray<T, DIM>& multiIndex) const
  {
    T rval = numerics::dot_product(multiIndex.begin(), m_strides.begin(), DIM);
    return rval;
  }

  //!@brief Convert flat index to multidimensional index.
  inline AXOM_HOST_DEVICE axom::StackArray<T, DIM> toMultiIndex(T flatIndex) const
  {
    axom::StackArray<T, DIM> multiIndex;
    for(int d = 0; d < DIM; ++d)
    {
      int dir = m_slowestDirs[d];
      multiIndex[dir] = flatIndex / m_strides[dir];
      flatIndex -= multiIndex[dir] * m_strides[dir];
    }
    return multiIndex;
  }
};

}  // end namespace axom

#endif  // QUEST_ARRAYINDEXER_HPP_
