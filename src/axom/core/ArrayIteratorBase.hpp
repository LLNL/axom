// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ARRAYITERATORBASE_HPP_
#define AXOM_ARRAYITERATORBASE_HPP_

#include "axom/core/IteratorBase.hpp"  // for Iterator

namespace axom
{
/// \name ArrayIteratorBase to iterate through Array-like types
/// @{

/**
 * \class   ArrayIteratorBase
 * \brief   An iterator type for Array-like types.
 *          Each increment operation advances the iterator to the next
 *          element in the Array-like.
 * \tparam ArrayType The type of the array to iterate over
 * \tparam ValueType The type of the array's elements
 */
template <typename ArrayType, typename ValueType>
class ArrayIteratorBase
  : public IteratorBase<ArrayIteratorBase<ArrayType, ValueType>, IndexType>
{
private:
  constexpr static bool ReturnConstRef = std::is_const<ValueType>::value;
  constexpr static bool FromArrayView =
    !std::is_const<typename ArrayType::RealConstT>::value;

public:
  // Iterator traits required to satisfy LegacyRandomAccessIterator concept
  // before C++20
  // See: https://en.cppreference.com/w/cpp/iterator/iterator_traits
  using difference_type = IndexType;
  using value_type = typename std::remove_cv<ValueType>::type;
  using reference = ValueType&;
  using pointer = ValueType*;
  using iterator_category = std::random_access_iterator_tag;

public:
  using ArrayPointerType =
    typename std::conditional<ReturnConstRef || FromArrayView,
                              const ArrayType*,
                              ArrayType*>::type;

  ArrayIteratorBase(IndexType pos, ArrayPointerType arr)
    : IteratorBase<ArrayIteratorBase<ArrayType, ValueType>, IndexType>(pos)
    , m_arrayPtr(arr)
  { }

  /**
   * \brief Returns the current iterator value
   */
  ValueType& operator*()
  {
    return (*m_arrayPtr)[IteratorBase<ArrayIteratorBase<ArrayType, ValueType>,
                                      IndexType>::m_pos];
  }

protected:
  /** Implementation of advance() as required by IteratorBase */
  void advance(IndexType n)
  {
    IteratorBase<ArrayIteratorBase<ArrayType, ValueType>, IndexType>::m_pos += n;
  }

protected:
  ArrayPointerType const m_arrayPtr;
};  // end of ArrayIteratorBase class

/// @}

} /* namespace axom */

#endif /* AXOM_ARRAYITERATORBASE_HPP_ */
