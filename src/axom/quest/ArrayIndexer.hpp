// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/StackArray.hpp"
#include "axom/core/numerics/matvecops.hpp"

namespace axom
{
/*!
  @brief Indexing into a multidimensional structured array.
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
        m_strides[d] = m_strides[d - 1] * lengths[d - 1];
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
        m_strides[d] = m_strides[d + 1] * lengths[d + 1];
      }
    }
  }

  //!@brief Constructor for arbitrary-stride indexing.
  ArrayIndexer(const axom::StackArray<T, DIM>& strides) : m_strides(strides)
  {
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
          std::swap(m_slowestDirs[s], m_slowestDirs[d]);
        }
      }
    }
  }

  //!@brief Index directions, ordered from slowest to fastest.
  inline AXOM_HOST_DEVICE axom::StackArray<std::uint16_t, DIM>& slowestDirs() const
  {
    return m_slowestDirs;
  }

  //!@brief Strides.
  inline AXOM_HOST_DEVICE axom::StackArray<axom::IndexType, DIM>& strides() const
  {
    return m_strides;
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
