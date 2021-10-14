// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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
 */
template <typename ArrayType>
class ArrayIteratorBase
  : public IteratorBase<ArrayIteratorBase<ArrayType>, IndexType>
{
public:
  ArrayIteratorBase(IndexType pos, ArrayType* arr)
    : IteratorBase<ArrayIteratorBase<ArrayType>, IndexType>(pos)
    , m_arrayPtr(arr)
  { }

  /**
   * \brief Returns the current iterator value
   */
  typename ArrayType::value_type& operator*()
  {
    return (
      *m_arrayPtr)[IteratorBase<ArrayIteratorBase<ArrayType>, IndexType>::m_pos];
  }

protected:
  /** Implementation of advance() as required by IteratorBase */
  void advance(IndexType n)
  {
    IteratorBase<ArrayIteratorBase<ArrayType>, IndexType>::m_pos += n;
  }

protected:
  ArrayType* const m_arrayPtr;
};  // end of ArrayIteratorBase class

/// @}

} /* namespace axom */

#endif /* AXOM_ARRAYITERATORBASE_HPP_ */
