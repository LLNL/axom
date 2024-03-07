// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_ARRAYINDEXER_HPP_
#define QUEST_ARRAYINDEXER_HPP_

#include "axom/slic.hpp"
#include "axom/core/StackArray.hpp"
#include "axom/core/numerics/matvecops.hpp"

#include <iostream>
#include <sstream>

namespace axom
{
/*!
  @brief Indicator for stride ordering.

  Multidimensional array data can be in row-major order, column-major order,
  or some arbitrarily permuted order.  Row and column major ordering are the
  same thing if the array is 1D.
*/
enum class ArrayStrideOrder : int
{
  ARBITRARY = 0,       // Neither row nor column
  ROW = 1,             // Row-major
  COLUMN = 2,          // Column-major
  BOTH = ROW | COLUMN  // 1D arrays are both
};

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
    @brief Constructor for row- or column-major indexing.
    @param [in] shape Shape of the array
    @param [in] arrayStrideOrder A order indicator from
                ArrayStrideOrder.
    @param [in] fastestStrideLength Stride in the fastest
                direction.
  */
  ArrayIndexer(const axom::StackArray<T, DIM>& shape,
               axom::ArrayStrideOrder arrayStrideOrder,
               int fastestStrideLength = 1)
  {
    initializeShape(shape, arrayStrideOrder, fastestStrideLength);
  }

  /*!
    @brief Constructor for a given order permutation.
    @param [in] shape Shape of the array
    @param [in] slowestDirs permutation vector, where
      slowestDirs[0] is the slowest direction and
      slowestDirs[DIM-1] is the fastest.
  */
  ArrayIndexer(const axom::StackArray<T, DIM>& shape,
               const axom::StackArray<std::uint16_t, DIM>& slowestDirs)
  {
    initializeShape(shape, slowestDirs);
  }

  /*!
    @brief Constructor for a given shape with the ordering of an
    existing ArrayIndexer.

    @param [in] shape Shape of the array
    @param [in] orderSource ArrayIndex to copy stride order
      from.
  */
  ArrayIndexer(const axom::StackArray<T, DIM>& shape,
               const axom::ArrayIndexer<T, DIM>& orderSource)
  {
    initializeShape(shape, orderSource);
  }

  /*!
    @brief Constructor for arbitrary-stride indexing.

    @param [i] strides Strides.  Must be unique when DIM > 1.
      If not unique, use default constructor and initializeStrides().

    @internal We could add the ArrayStrideOrder preference to this constructor
    to handle the degenerate case of non-unique strides.  But that would
    clash with the more prevalent usage of constructing from the array's
    shape.
  */
  ArrayIndexer(const axom::StackArray<T, DIM>& strides) : m_strides(strides)
  {
    initializeStrides(strides);
  }

  /*!
    @brief Default constructor

    Object must be initialized before use.
  */
  ArrayIndexer() = default;

  /*!
    @brief Initialize for row- or column-major indexing.
    @param [in] shape Shape of the array
    @param [in] arrayStrideOrder An order indicator from
                ArrayStrideOrder.
    @param [in] fastestStrideLength Stride in the fastest
                direction.
  */
  inline AXOM_HOST_DEVICE void initializeShape(const axom::StackArray<T, DIM>& shape,
                                               ArrayStrideOrder arrayStrideOrder,
                                               int fastestStrideLength = 1)
  {
    SLIC_ASSERT(arrayStrideOrder == ArrayStrideOrder::COLUMN ||
                arrayStrideOrder == ArrayStrideOrder::ROW ||
                (DIM == 1 && arrayStrideOrder == ArrayStrideOrder::BOTH));
    if(arrayStrideOrder == ArrayStrideOrder::ROW)
    {
      for(int d = 0; d < DIM; ++d)
      {
        m_slowestDirs[d] = DIM - 1 - d;
      }
      m_strides[0] = fastestStrideLength;
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
      m_strides[DIM - 1] = fastestStrideLength;
      for(int d = DIM - 2; d >= 0; --d)
      {
        m_strides[d] = m_strides[d + 1] * shape[d + 1];
      }
    }
  }

  /*!
    @brief Initialize for a given order permutation.
    @param [in] shape Shape of the array
    @param [in] slowestDirs permutation vector, where
      slowestDirs[0] is the slowest direction and
      slowestDirs[DIM-1] is the fastest.
  */
  inline AXOM_HOST_DEVICE void initializeShape(
    const axom::StackArray<T, DIM>& shape,
    const axom::StackArray<std::uint16_t, DIM>& slowestDirs)
  {
    SLIC_ASSERT(isPermutation(slowestDirs));
    m_slowestDirs = slowestDirs;
    m_strides[m_slowestDirs[DIM - 1]] = 1;
    for(int d = DIM - 2; d >= 0; --d)
    {
      int dir = m_slowestDirs[d];
      int fasterDir = m_slowestDirs[d + 1];
      m_strides[dir] = m_strides[fasterDir] * shape[fasterDir];
    }
  }

  /*!
    @brief Initialize for a given shape with the ordering of an
    existing ArrayIndexer.

    @param [in] shape Shape of the array
    @param [in] orderSource ArrayIndex to copy stride order
      from.
  */
  inline AXOM_HOST_DEVICE void initializeShape(
    const axom::StackArray<T, DIM>& shape,
    const axom::ArrayIndexer<T, DIM>& orderSource)
  {
    initializeShape(shape, orderSource.slowestDirs());
  }

  /*!
    @brief Initialize for arbitrary-stride indexing.

    @param [i] strides Strides.  Must be unique when DIM > 1.
      If not satisfied, you must use one of the other initializers.
  */
  inline AXOM_HOST_DEVICE void initializeStrides(
    const axom::StackArray<T, DIM>& strides)
  {
    if(DIM > 1 && !stridesAreUnique(strides))
    {
#if !defined(AXOM_DEVICE_CODE)
      std::ostringstream os;
      os << "(";
      for(int d = 0; d < DIM - 1; ++d)
      {
        os << strides[d] << ",";
      }
      os << strides[DIM - 1] << ")";
      std::cerr << "ERROR: ArrayIndexer: Non-unique strides " << os.str() << ".\n"
                << "Likely, multi-dim array shape is 1 in some direction.\n"
                << "Impossible to compute index ordering.\n"
                << "Please use a different ArrayIndexer initializer.\n";
#endif
      utilities::processAbort();
    }

    // 2nd argument doesn't matter because strides are unique.
    initializeStrides(strides, axom::ArrayStrideOrder::COLUMN);
  }

  /*!
    @brief Initialize for arbitrary-stride indexing,
    with ordering preference for non-unique strides.

    @param [i] strides Strides.
    @param [i] orderPref Ordering preference value
      (from ArrayStrideOrder) if strides are non-unique.
  */
  inline AXOM_HOST_DEVICE void initializeStrides(
    const axom::StackArray<T, DIM>& strides,
    ArrayStrideOrder orderPref)
  {
    SLIC_ASSERT(orderPref == axom::ArrayStrideOrder::COLUMN ||
                orderPref == axom::ArrayStrideOrder::ROW);

    m_strides = strides;
    for(int d = 0; d < DIM; ++d)
    {
      m_slowestDirs[d] =
        orderPref == axom::ArrayStrideOrder::COLUMN ? d : DIM - 1 - d;
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

  static AXOM_HOST_DEVICE bool stridesAreUnique(const axom::StackArray<T, DIM>& strides)
  {
    bool repeats = false;
    for(int d = 1; d < DIM; ++d)
    {
      for(int e = 0; e < d; ++e)
      {
        repeats |= strides[d] == strides[e];
      }
    }
    return !repeats;
  }

  inline AXOM_HOST_DEVICE bool operator==(const ArrayIndexer& other) const
  {
    return m_slowestDirs == other.m_slowestDirs && m_strides == other.m_strides;
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

  //!@brief Whether a StackArray represents a permutation.
  bool isPermutation(const axom::StackArray<std::uint16_t, DIM>& v)
  {
    // v is a permutation if all its values are unique and in [0, DIM).
    axom::StackArray<bool, DIM> found;
    for(int d = 0; d < DIM; ++d)
    {
      found[d] = false;
    }
    for(int d = 0; d < DIM; ++d)
    {
      if(v[d] < 0 || v[d] >= DIM)
      {
        return false;  // Out of range.
      }
      if(found[v[d]] == true)
      {
        return false;  // Repeated index.
      }
      found[v[d]] = true;
    }
    return true;
  }

  /*!
    @brief Get the stride order (row- or column-major, or something else).

    @return Value from ArrayStrideOrder, indicating column order,
       row order, both column and row (1D only) or arbitrary order.
  */
  inline AXOM_HOST_DEVICE ArrayStrideOrder getStrideOrder() const
  {
    int ord = int(ArrayStrideOrder::BOTH);
    for(int d = 0; d < DIM - 1; ++d)
    {
      ord &= m_slowestDirs[d] < m_slowestDirs[d + 1]
        ? int(ArrayStrideOrder::COLUMN)
        : int(ArrayStrideOrder::ROW);
    }
    static ArrayStrideOrder s_intToOrder[4] = {ArrayStrideOrder::ARBITRARY,
                                               ArrayStrideOrder::ROW,
                                               ArrayStrideOrder::COLUMN,
                                               ArrayStrideOrder::BOTH};
    return s_intToOrder[ord];
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

  friend std::ostream& operator<<(std::ostream& os, const ArrayIndexer& a)
  {
    os << "ArrayIndexer: strides=(" << a.m_strides[0];
    for(int d = 1; d < DIM; ++d)
    {
      os << ',' << a.m_strides[d];
    }
    os << ") slowestDirs=(" << a.m_slowestDirs[0];
    for(int d = 1; d < DIM; ++d)
    {
      os << ',' << a.m_slowestDirs[d];
    }
    os << ')';
    return os;
  }
};

}  // end namespace axom

#endif  // QUEST_ARRAYINDEXER_HPP_
